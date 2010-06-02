###
# Pod Documentation
###

=head1 NAME

Breseq::Fastq.pm

=head1 SYNOPSIS

Module for reading and writing fastq files more rapidly than BioPerl.

=head1 DESCRIPTION

=head1 TODO

Read translation table from GenBank file

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

use Bio::Seq::RichSeq;
use Breseq::Settings;
use Breseq::Shared;
use Data::Dumper;


##this is a version of the next subroutine that only loads the features for annotation
sub load_ref_seq_info
{
	##summaryand create_fasta are optional
	my ($settings, $summary, $create_fasta) = @_;
			
	print STDERR "Loading reference sequences...\n";

	$summary->{sequence_conversion}->{total_reference_sequence_length} = 0 if (defined $summary);
	
	my $ref_seq_info;
	my $i = 0;
	my $s;

	# get two pieces of information from $settings
	my @genbank_file_names = $settings->file_name('reference_genbank_file_names'); 
	my @junction_only_genbank_file_names = $settings->file_name('junction_only_reference_genbank_file_names'); 

	my %junction_only_hash;
	foreach my $jo (@junction_only_genbank_file_names)
	{
		$junction_only_hash{$jo} = 1;
	}
	
	## optionally create FASTA for writing to...
	my $fasta_o;
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	if ($create_fasta)
	{
		$fasta_o = Bio::SeqIO->new( -format => "FASTA", -file => ">$reference_fasta_file_name");
	}	

	foreach my $genbank_file_name (@genbank_file_names, @junction_only_genbank_file_names)
	{
		## open this GenBank file
		my $ref_i = Bio::SeqIO->new( -file => "<$genbank_file_name");
		print STDERR "  Loading File::$genbank_file_name\n";

		my $junction_only = $junction_only_hash{$genbank_file_name};

		while (my $ref_seq = $ref_i->next_seq)
		{			
			my $seq_id = $ref_seq->id;
			print STDERR "    Sequence::$seq_id loaded.\n";
			
			$s->{$seq_id}->{seq_id} = $seq_id;
			$s->{$seq_id}->{length} = $ref_seq->length;
			$s->{$seq_id}->{definition} = $ref_seq->desc;
			$s->{$seq_id}->{version} = $ref_seq->display_id;			
			
			$summary->{sequence_conversion}->{total_reference_sequence_length} += $ref_seq->length if (defined $summary);
			
			## it would be nice to get rid of storing the whole genome in memory
			$ref_seq_info->{bioperl_ref_seqs}->{$seq_id} = $ref_seq;

			##it is much faster to use substr to create the lists for nt comparisons than BioPerl trunc
			
			#load the reference sequence, uppercase it, and scrub degenerate bases to N's 
			$ref_seq_info->{ref_strings}->{$seq_id} = $ref_seq->seq;
			$ref_seq_info->{ref_strings}->{$seq_id} = uc($ref_seq_info->{ref_strings}->{$seq_id});
			$ref_seq_info->{ref_strings}->{$seq_id} =~ s/[^ATCG]/N/g;
			$ref_seq->seq($ref_seq_info->{ref_strings}->{$seq_id});						
			$ref_seq_info->{seq_order}->{$seq_id} = $i++;
			
			## optionally write out FASTA -- AFTER scrubbing sequence
			if ($create_fasta)
			{
				$fasta_o->write_seq($ref_seq);
			}
			
			if (!$junction_only)
			{
				push @{$ref_seq_info->{seq_ids}}, $seq_id;
			}
			else
			{
				push @{$ref_seq_info->{junction_only_seq_ids}}, $seq_id;
			}

			##load the genbank record
			my @Feature_List = $ref_seq->get_SeqFeatures();
			my @gene_list;
			my @repeat_list;
			
			FEATURE: foreach my $Feature ( @Feature_List ) 
			{ 	
				if ($Feature->primary_tag eq 'repeat_region')
				{	
					my $r;
					$r->{name} = get_tag($Feature, "mobile_element");
					$r->{name} =~ s/insertion sequence:// if ($r->{name});
					
					## don't add unnamed ones to list...
					#$r->{name} = "unknown_repeat" if (!$r->{name});
					next FEATURE if (!$r->{name});
					
					$r->{product} = "repeat region";
					$r->{start} = $Feature->start;
					$r->{end} = $Feature->end;
					$r->{strand} = $Feature->strand;
					push @repeat_list, $r;
					next FEATURE;
				}
				
				## add additional information to the last 
				if (   ($Feature->primary_tag eq 'CDS') 
					|| ($Feature->primary_tag eq 'tRNA') 
					|| ($Feature->primary_tag eq 'rRNA') )

				{
					#Add information to last gene
					my $gene;
					$gene->{name} = get_tag($Feature, "gene");
					$gene->{name} = get_tag($Feature, "locus_tag") if (!$gene->{name});
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
			
			$ref_seq_info->{gene_lists}->{$seq_id} = \@gene_list;
			$ref_seq_info->{repeat_lists}->{$seq_id} = \@repeat_list;	
		}
	}
			
	#save statistics if requested
	$summary->{sequence_conversion}->{reference_sequences} = $s if (defined $summary);

	## create SAM faidx
	Breseq::Shared::system("samtools faidx $reference_fasta_file_name", 1) if ($create_fasta);
			
	return $ref_seq_info;
}

sub annotate_1_mutation
{
	my ($ref_seq_info, $mut, $start, $end, $repeat_override) = @_;

	## this could be moved to the object
	my $intergenic_seperator = "/";

	## initialize everything, even though we don't always use it
	$mut->{aa_position} = "";
	$mut->{aa_ref_seq} = "";
	$mut->{aa_new_seq} = "";
	$mut->{codon_position} = "";
	$mut->{codon_ref_seq} = "";
	$mut->{codon_new_seq} = "";		
	$mut->{gene_name} = "";
	$mut->{gene_position} = "";
	$mut->{gene_product} = "";

	my $seq_id = $mut->{seq_id};
	my $gene_list_ref = $ref_seq_info->{gene_lists}->{$seq_id};
	my $repeat_list_ref = $ref_seq_info->{repeat_lists}->{$seq_id};
	my $ref_seq = $ref_seq_info->{bioperl_ref_seqs}->{$seq_id};		

	my $size = $end - $start + 1;

	my ($prev_gene, $next_gene) = (undef, undef);
	my @within_genes = ();
	my @between_genes = ();
	my @inside_left_genes = ();
	my @inside_right_genes = ();
		
	my $repeat_region;
	if ($repeat_override)
	{
		die if ($start != $end);
		$repeat_region = get_overlapping_feature($repeat_list_ref, $start);
		push @within_genes, $repeat_region if ($repeat_region);
	}
	
	if (!$repeat_region)
	{	
		my ($within_genes_list_ref, $between_genes_list_ref, $inside_left_genes_list_ref, $inside_right_genes_list_ref);
		($prev_gene, $next_gene, $within_genes_list_ref, $between_genes_list_ref, $inside_left_genes_list_ref, $inside_right_genes_list_ref) 
			= find_nearby_genes($gene_list_ref, $start, $end);	
			
		@within_genes = @$within_genes_list_ref;
		@between_genes = @$between_genes_list_ref;
		@inside_left_genes = @$inside_left_genes_list_ref;
		@inside_right_genes = @$inside_right_genes_list_ref;
	}


	## Mutation is noncoding
	if (scalar(@within_genes) + scalar(@between_genes) == 0)
	{			
		$mut->{snp_type} = "IG";

		$mut->{gene_name} .= (defined $prev_gene) ? $prev_gene->{name} : "&minus;";
		$mut->{gene_name} .= $intergenic_seperator;
		$mut->{gene_name} .= (defined $next_gene) ? $next_gene->{name} : "&minus;";

		if (defined $prev_gene)
		{
			$mut->{gene_position} .= "intergenic&nbsp;(";
			$mut->{gene_position} .= ($prev_gene->{strand} == +1) ? "+" : "-";
			$mut->{gene_position} .= $start - $prev_gene->{end};
		}
		else
		{
			$mut->{gene_position} .= "intergenic&nbsp;(&minus;";
		}
		$mut->{gene_position} .= $intergenic_seperator;
		if (defined $next_gene)
		{
			$mut->{gene_position} .= ($next_gene->{strand} == +1) ? "-" : "+";
			$mut->{gene_position} .= $next_gene->{start} - $end;
		}
		else
		{
			$mut->{gene_position} .= "&minus;";
		}
		$mut->{gene_position} .= ")";

		$mut->{gene_product} .= (defined $prev_gene) ? $prev_gene->{product} : "&minus;";
		$mut->{gene_product} .= $intergenic_seperator;			
		$mut->{gene_product} .= (defined $next_gene) ? $next_gene->{product} : "&minus;";

		return $mut;
	}
	## Mutation is completely within genes
	elsif (scalar @within_genes > 0)
	{
### TO DO: It can be within multiple genes, in which case we need to annotate
### the change it causes in each reading frame UGH! YUCKY!				
### FOR NOW: just take the first of the within genes...				
		my $gene = $within_genes[0];
		$mut->{gene_name} = $gene->{name};
		$mut->{gene_product} = $gene->{product};

		my $within_gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};	

		if ($start == $end)
		{
			$mut->{gene_position} = abs($start-$within_gene_start) + 1; 
		}
		else
		{
			my $gene_start = abs($start-$within_gene_start) + 1; 
			my $gene_end = abs($end-$within_gene_start) + 1; 
			$mut->{gene_position} = ($gene_start < $gene_end) ? "$gene_start&minus;$gene_end" : "$gene_end&minus;$gene_start";
		}

		my $gene_nt_size = $gene->{end} - $gene->{start} + 1;

		## ...but the gene is a pseudogene or an RNA gene
		if ( (!$gene->{translation} || $gene->{pseudogene}) )
		{
			$mut->{snp_type} = "NC";
			$mut->{gene_position} = "noncoding ($mut->{gene_position}/$gene_nt_size&nbsp;nt)";
			return $mut;
		}	

		if ($mut->{type} ne 'SNP')
		{					
			$mut->{gene_position} = "coding&nbsp;($mut->{gene_position}/$gene_nt_size&nbsp;nt)";
			return $mut;
		}

	    ## determine the old and new translation of this codon  
   		$mut->{aa_position} = int(($mut->{gene_position}-1)/3) + 1;
		$mut->{codon_position} = abs($start-$within_gene_start) % 3;

		my $codon_seq = ($gene->{strand} == +1) ?
			$ref_seq->trunc($gene->{start} + 3 * ($mut->{aa_position}-1),$gene->{start} + 3 * $mut->{aa_position} - 1) :
			$ref_seq->trunc($gene->{end} - 3 * $mut->{aa_position}+1,$gene->{end} - 3 * ($mut->{aa_position}-1))->revcom;

		$mut->{codon_ref_seq} = $codon_seq->seq();
		$mut->{aa_ref_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();

		$mut->{codon_new_seq} = $codon_seq->seq();
		#remember to revcom the change if gene is on opposite strand
		substr($mut->{codon_new_seq}, $mut->{codon_position}, 1) = ($gene->{strand} == +1) ? $mut->{new_seq} : revcom($mut->{new_seq});
		$codon_seq->seq($mut->{codon_new_seq});
		$mut->{aa_new_seq} = $codon_seq->translate( undef, undef, undef, 11 )->seq();

	    $mut->{snp_type} = ($mut->{aa_ref_seq} ne $mut->{aa_new_seq}) ? "NS" : "S";
	}

	##The mutation actually contains several genes
	elsif (scalar @between_genes > 0)
	{
		my @gene_list = ( map({ "<i>[" . $_->{name} . "]</i>" } @inside_left_genes),
						  map({ "<i>" . $_->{name} . "</i>" } @between_genes),
						  map({ "<i>[" . $_->{name} ."]</i>" } @inside_right_genes) );

		$mut->{gene_product} = join ("&minus", @gene_list);


		if (scalar @gene_list == 1)
		{
			$mut->{gene_name} = $gene_list[0];
		}
		else
		{
			$mut->{gene_name} = $gene_list[0] . "&minus" . $gene_list[-1];
		}
	}
	
	return $mut;
}

sub annotate_mutations
{
	my ($ref_seq_info, $gd) = @_;

	##keep track of other mutations that affect SNPs
	##because we may double-hit a codon

#TODO: the proper way to do this is to create list of SNPs that have been hit
# hashed by gene protein accession ID and AA position within gene
# and have the annotation point to them (and back at them)
# so that the codon will be correctly updated with all changes and we can notify the
# changes that their SNP_type is not really SNP, but multiple hit SNP.
	
	my $snp_hits_hash;
	
	MUT: foreach my $mut ($gd->list)
	{		
		
		if ($mut->{type} eq 'SNP')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
		}
		elsif ($mut->{type} eq 'SUB')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + length($mut->{ref_seq}) - 1);
		}
		elsif ($mut->{type} eq 'DEL')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position} + $mut->{size} - 1);
		}
		elsif ($mut->{type} eq 'INS')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
		}
		elsif ($mut->{type} eq 'MOB')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{start}, $mut->{end});
		}
		elsif ($mut->{type} eq 'JC')
		{
			annotate_1_mutation($ref_seq_info, $mut->{_side_1}, $mut->{side_1_position}, $mut->{side_1_position}, 1);
			annotate_1_mutation($ref_seq_info, $mut->{_side_2}, $mut->{side_2_position}, $mut->{side_2_position}, 1);
			
			print Dumper($mut);
		}
		elsif ($mut->{type} eq 'RA')
		{
			annotate_1_mutation($ref_seq_info, $mut, $mut->{position}, $mut->{position});
		}
	}
}

sub revcom
{
	my ($seq) = @_;
	$seq =~ tr/ATCG/TAGC/;
	return reverse $seq;
}

sub annotate_deletions
{
	my ($settings, $summary, $ref_seq_info, $deletions_ref) = @_;	
	
	foreach my $deletion (@$deletions_ref)
	{
		#set up to use correct reference sequence
		my $seq_id = $deletion->{seq_id};
		my $gene_list_ref = $ref_seq_info->{gene_lists}->{$seq_id};
		my $repeat_list_ref = $ref_seq_info->{repeat_lists}->{$seq_id};
		
		my $gene_string = '';
		my $feature_start = get_overlapping_feature($repeat_list_ref, $deletion->{start});
		$feature_start = get_overlapping_feature($gene_list_ref, $deletion->{start}) if (!$feature_start);
	
		my $feature_end = get_overlapping_feature($repeat_list_ref, $deletion->{end}, );
		$feature_end = get_overlapping_feature($gene_list_ref, $deletion->{end}, ) if (!$feature_end);
	
		if ($feature_start)
		{
			$gene_string .= "[$feature_start->{name}]";
		}
	
		my $start_bound = ($feature_start) ? $feature_start->{end}-25 : $deletion->{start};
		my $end_bound = ($feature_end) ? $feature_end->{start} : $deletion->{end};

		foreach my $gene (@$gene_list_ref)
		{
			last if ($gene->{start} >= $end_bound);
		
			if ($gene->{start} > $start_bound)
			{
				$gene_string .= " " if ($gene_string);
				$gene_string .= "$gene->{name}";
			}
		}
	
		## only add the end gene if it isn't the start gene
		if ($feature_end && (!$feature_start || ($feature_start != $feature_end)))
		{
			$gene_string .= " " if ($gene_string);
			$gene_string .= "[$feature_end->{name}]";
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
		## This additional information is used for the complex reference line.
		## Note that it is completely determined by the original candidate junction sequence 
		## positions and overlap: alignment_pos and alignment_overlap.
		my $alignment_reference_info_1 = { 
			truncate_end => $item->{flanking_left}, 
			ghost_end => $item->{side_1}->{alignment_pos}, 
			ghost_strand => $item->{side_1}->{strand},
			ghost_seq_id => $item->{side_1}->{seq_id}
		};
		$alignment_reference_info_1->{truncate_end} += $item->{alignment_overlap} if ($item->{alignment_overlap} > 0);

		my $alignment_reference_info_2 = { 
			truncate_start => $item->{flanking_left}+1, 
			ghost_start => $item->{side_2}->{alignment_pos}, 
			ghost_strand => $item->{side_2}->{strand},
			ghost_seq_id => $item->{side_2}->{seq_id}
		};
		$alignment_reference_info_2->{truncate_start} -= $item->{alignment_overlap} if ($item->{alignment_overlap} < 0);
		
		push @{$item->{alignment_reference_info_list}}, $alignment_reference_info_1, $alignment_reference_info_2;
		$item->{alignment_empty_change_line} = 1;

		#add gene information for each end
		$item->{hybrid} = $item;
		foreach my $key ('side_1', 'side_2')
		{			
			##create circular reference to self so we can print table at the top of the alignment
			$item->{$key}->{hybrid} = $item;
			
			my ($prev_gene, $gene, $next_gene) 
				= Breseq::ReferenceSequence::find_nearby_genes($ref_seq_info->{gene_lists}->{$item->{$key}->{seq_id}}, $item->{$key}->{start});		

			## noncoding
			if (!$gene)
			{
				$item->{$key}->{name}->{name} .= ($prev_gene && $prev_gene->{name}) ? $prev_gene->{name} : '-';
				$item->{$key}->{name}->{name} .= "/";
				$item->{$key}->{name}->{name} .= ($next_gene && $next_gene->{name}) ? $next_gene->{name} : '-';
	
				if (defined $prev_gene)
				{
					$item->{$key}->{name}->{position} .= ($prev_gene->{strand} == +1) ? "+" : "-";
					$item->{$key}->{name}->{position} .= ($item->{$key}->{start} - $prev_gene->{end});
				}
				$item->{$key}->{name}->{gene_position} .= "/";
				if (defined $next_gene)
				{
					$item->{$key}->{name}->{position} .= ($next_gene->{strand} == +1) ? "-" : "+";
					$item->{$key}->{name}->{position} .= ($next_gene->{start} - $item->{$key}->{end});
				}
		
				$item->{$key}->{name}->{product} .= ($prev_gene && $prev_gene->{product}) ? $prev_gene->{product} : '-';
				$item->{$key}->{name}->{product} .= "/";			
				$item->{$key}->{name}->{product} .= ($next_gene && $next_gene->{product}) ? $next_gene->{product} : '-';

				$item->{$key}->{name}->{interval} .= $prev_gene->{end}+1 if $prev_gene;
				$item->{$key}->{name}->{interval} .= "/"; 
				$item->{$key}->{name}->{interval} .= $next_gene->{start}-1 if $next_gene; 
			}
			## coding
			else
			{
				$item->{$key}->{name}->{name} = $gene->{name};
				$item->{$key}->{name}->{product} = $gene->{product};
				my $gene_start = ($gene->{strand} == +1) ? $gene->{start} : $gene->{end};	
				$item->{$key}->{name}->{position} = abs($item->{$key}->{start}-$gene_start) + 1;
				$item->{$key}->{name}->{interval} = ($gene->{strand} == +1) ? "$gene->{start}-$gene->{end}" : "$gene->{end}-$gene->{start}"; 
			}
		
			## Determine IS elements
			## Is it within an IS or near the boundary of an IS in the direction leading up to the junction?			
			if (my $is = Breseq::ReferenceSequence::find_closest_repeat_region($item->{$key}, $ref_seq_info->{repeat_lists}->{$item->{$key}->{seq_id}}, 200, $item->{$key}->{strand}))
			{
				$item->{$key}->{is}->{name} = $is->{name};
				$item->{$key}->{is}->{interval} = ($is->{strand} == +1) ? "$is->{start}-$is->{end}" : "$is->{end}-$is->{start}"; 
				$item->{$key}->{is}->{product} = $is->{product};
			}
			$item->{$key}->{annotate_key} = (defined $item->{$key}->{is}) ? 'is' : 'gene';			
		}
		print Dumper($item) if ($verbose);
	}

	if (!$settings->{sort_junctions_by_score})
	{
		sub by_hybrid
		{
			my $a_pos = (defined $a->{side_1}->{is}) ? $a->{side_2}->{start} : $a->{side_1}->{start};
			my $b_pos = (defined $b->{side_1}->{is}) ? $b->{side_2}->{start} : $b->{side_1}->{start};
	
			my $a_seq_order = (defined $a->{side_1}->{is}) ? $ref_seq_info->{seq_order}->{$a->{side_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$a->{side_1}->{seq_id}};
			my $b_seq_order = (defined $b->{side_1}->{is}) ? $ref_seq_info->{seq_order}->{$b->{side_2}->{seq_id}} : $ref_seq_info->{seq_order}->{$b->{side_1}->{seq_id}};		
	
			return (($a_seq_order <=> $b_seq_order) || ($a_pos <=> $b_pos));
		}
		@$rearrangements_ref = sort by_hybrid @$rearrangements_ref;
	}

	
	###
	# IS insertion overlap correction
	#
	# For these the coordinates may have been offset incorrectly initially (because both sides of the junction may look unique)
	# The goal is to offset through positive overlap to get as close as possible to the ends of the IS
	###
	foreach my $j (@$rearrangements_ref)
	{	
		
		sub add_is_coords_from_interval
		{
			my ($c) = @_;
			return if (!defined $c->{is}); 
			
			my ($is_start, $is_end) = split /-/, $c->{is}->{interval};
			$c->{is}->{strand} = ($is_start < $is_end) ? +1 : -1; 
			$c->{is}->{start} = ($is_start < $is_end) ? $is_start : $is_end; 
			$c->{is}->{end} = ($is_start < $is_end) ? $is_end : $is_start;
		}
		
		add_is_coords_from_interval($j->{side_1});
		add_is_coords_from_interval($j->{side_2});
		
		$j->{side_1}->{read_side} = -1;
		$j->{side_2}->{read_side} = +1;
				
		## Determine which side of the junction is the IS and which is unique
		## these point to the correct initial interval...
		if (defined $j->{side_1}->{is})
		{
			if (abs($j->{side_1}->{is}->{start} - $j->{side_1}->{start}) <= 20)
			{
				$j->{is_interval} = $j->{side_1};
				$j->{is_interval}->{is}->{side_key} = 'start';
			}
			elsif (abs($j->{side_1}->{is}->{end} - $j->{side_1}->{start}) <= 20 )
			{
				$j->{is_interval} = $j->{side_1};
				$j->{is_interval}->{is}->{side_key} = 'end';
			}
			$j->{unique_interval} = $j->{side_2};
		}
		
		if (!defined $j->{is} && defined $j->{side_2}->{is})
		{
			if (abs($j->{side_2}->{is}->{start} - $j->{side_2}->{start}) <= 20)
			{
				$j->{is_interval} = $j->{side_2};
				$j->{is_interval}->{is}->{side_key} = 'start';
			}
			elsif (abs($j->{side_2}->{is}->{end} - $j->{side_2}->{start}) <= 20 )
			{
				$j->{is_interval} = $j->{side_2};
				$j->{is_interval}->{is}->{side_key} = 'end';
			}
			$j->{unique_interval} = $j->{side_1};
		}
		
		## Ah, we don't have an IS, we are done
		next if (!defined $j->{is_interval});
		
		## Ah, there is no overlap to play with, we are done
		next if ($j->{overlap} <= 0);
		
		## The following code implies $j->{overlap} > 0
				
		### first, adjust the repetitive sequence boundary to get as close to the IS as possible
		my $move_dist = abs($j->{is_interval}->{start} - $j->{is_interval}->{is}->{$j->{is_interval}->{is}->{side_key}});
		$move_dist = $j->{overlap} if ($move_dist > $j->{overlap});
		$j->{is_interval}->{start} += $j->{is_interval}->{strand} * $move_dist;
		$j->{overlap} -= $move_dist;
		$j->{is_interval}->{end} = $j->{is_interval}->{start};
		
		### second, adjust the unique sequence side with any remaining overlap
		$j->{unique_interval}->{start} += $j->{unique_interval}->{strand} * $j->{overlap};	
		$j->{unique_interval}->{end} = $j->{unique_interval}->{start};
				
		$j->{overlap} = 0;
	}
		
}


sub get_overlapping_feature
{
	my ($feature_list_ref, $pos) = @_;
	foreach my $f (@$feature_list_ref)
	{
		return $f if ( ($pos >= $f->{start}) && ($pos <= $f->{end}) );
	}
	return undef;
}

sub find_nearby_genes
{
	my ($gene_list_ref, $pos_1, $pos_2) = @_;
	$pos_2 = $pos_1 if (!defined $pos_2);

#	print "$pos_1, $pos_2\n";

	my (@within_genes, @between_genes, @inside_left_genes, @inside_right_genes, $prev_gene, $next_gene);
	GENE: for (my $i=0; $i < scalar @$gene_list_ref; $i++)
	{
		my $test_gene = $gene_list_ref->[$i];

		if ($test_gene->{end} < $pos_1)
		{
			$prev_gene = $test_gene;
		}
		
		if (  ($test_gene->{start} <= $pos_1) && ($test_gene->{end} >= $pos_1) 
		   && ($test_gene->{start} <= $pos_2) && ($test_gene->{end} >= $pos_2) )
		{
			push @within_genes, $test_gene;
#			print "^ $test_gene->{name}\n";			
		}
		elsif ( ($test_gene->{start} <= $pos_1) && ($test_gene->{end} >= $pos_1) )
		{
			push @inside_left_genes, $test_gene;
		}
		elsif ( ($test_gene->{start} >= $pos_2) && ($test_gene->{end} <= $pos_2) )
		{
			push @inside_right_genes, $test_gene;
		}
		elsif ( ($test_gene->{start} >= $pos_1) && ($test_gene->{end} <= $pos_2) )
		{	
			push @between_genes, $test_gene;
#			print ">< $test_gene->{name}\n";
		}
		#We've passed the changes, so it is in the previous intergenic space
		if ($test_gene->{start} > $pos_2)
		{
			$next_gene = $test_gene;
			last GENE;
		}
	}

#	print "$prev_gene->{name} || $next_gene->{name}\n";

	return ($prev_gene, $next_gene, \@within_genes, \@between_genes, \@inside_left_genes, \@inside_right_genes);
}

sub find_closest_repeat_region
{
	my ($position, $repeat_list_ref, $max_distance, $direction) = @_;

	my $is = undef;
	my $best_distance;
	IS: for (my $i=0; $i < scalar @$repeat_list_ref; $i++)
	{
		my $test_is = $repeat_list_ref->[$i];
				
		#count within the IS element as zero distance
		#if this happens then we are immediately done
		if ( ($test_is->{start} <= $position) && ($test_is->{end} >= $position) )
		{	
			return $test_is;
		}
		
		#otherwise calculate the distance
		#keep if less than max_distance, in the correct direction, and the best found so far 
		my $test_distance = ($direction == -1) ? $position - $test_is->{end} : $test_is->{start} - $position;
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