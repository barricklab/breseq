###
# Pod Documentation
###

=head1 NAME
Breseq::AlignmentOutput.pm

=head1 SYNOPSIS

Module for making alignments for display from SAM databases.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jeffrey@barricklab.org>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

use strict;

package Breseq::AlignmentOutput;
use vars qw(@ISA);
use Bio::Root::Root;
@ISA = qw( Bio::Root::RootI );

use List::Util qw(shuffle);
use CGI qw/:standard *table *Tr *code *td start_b end_b start_i end_i/;

use Bio::DB::Sam;
use Breseq::Fastq;
use Breseq::Shared;
use Data::Dumper;

## "rgb(0,255,255)" = cyan

our $base_colors_hash = {
	'G' => [ "rgb(255,255,0)", "rgb(230,230,230)", "rgb(210,210,210)", "rgb(140,140,140)", "rgb(70,70,70)", "rgb(0,0,0)" ],
	'C' => [ "rgb(255,255,0)", "rgb(160,160,255)", "rgb(120,120,255)", "rgb(60,60,255)", "rgb(0,0,255)", "rgb(0,0,150)" ],
	'A' => [ "rgb(255,255,0)", "rgb(255,210,210)", "rgb(255,180,180)", "rgb(255,100,100)", "rgb(255,20,20)", "rgb(200,0,0)" ],
	'T' => [ "rgb(255,255,0)", "rgb(210,255,210)", "rgb(180,255,180)", "rgb(100,255,100)", "rgb(20,255,20)", "rgb(0,200,0)" ],
	'N' => [ "rgb(128,0,128)", "rgb(128,0,128)", "rgb(128,0,128)", "rgb(128,0,128)", "rgb(128,0,128)", "rgb(128,0,128)" ],
};

sub new
{
	my($caller,@args) = @_;
	my $class = ref($caller) || $caller;
	my $self = new Bio::Root::Root($caller, @args);

	bless ($self, $class);
	($self->{maximum_to_align}) = $self->Bio::Root::RootI::_rearrange([qw(MAXIMUM_TO_ALIGN)], @args);
	$self->{maximum_to_align} = 200 if (!defined $self->{maximum_to_align});
	$self->{maximum_to_make_alignment} = 100000 if (!defined $self->{maximum_to_make_alignment});

	$self->{header_style_string} = '';
	$self->{header_style_string} .= "\.NC {color: rgb(0,0,0); background-color: rgb(255,255,255)}\n"; #no color
	$self->{header_style_string} .= "\.UN {color: rgb(120,120,120); background-color: rgb(255,255,255)}\n"; #unaligned
	foreach my $key (keys %$base_colors_hash)
	{
		for (my $i=0; $i<scalar @{$base_colors_hash->{$key}}; $i++)
		{
			
			if ($i>0)
			{
				$self->{header_style_string} .= "\.$key$i \{color: rgb(255,255,255); background-color: $base_colors_hash->{$key}->[$i]\}\n";
			}
			else
			{
				$self->{header_style_string} .= "\.$key$i \{color: rgb(120,120,120); background-color: $base_colors_hash->{$key}->[$i]\}\n";
			}
		}
	}
	
	$self->{no_color_index} = scalar(@{$base_colors_hash->{'G'}}) - 1;

	return $self; 
}


sub text_alignment
{
	my ($self, $bam_path, $fasta_path, $region, $options) = @_;
	$options->{text} = 1;
	my $alignment_info = $self->create_alignment($bam_path, $fasta_path, $region, $options);
	return 'No reads uniquely aligned to position.' if (!defined $alignment_info);
	
	my $aligned_reads = $alignment_info->{aligned_reads};
	my @aligned_references = @{$alignment_info->{aligned_references}};
	
	my $output = '';
			
	my @sorted_keys = sort { -($aligned_reads->{$a}->{aligned_bases} cmp $aligned_reads->{$b}->{aligned_bases}) } keys %$aligned_reads;
	
	foreach my $aligned_reference (@aligned_references)
	{
		$output .= _text_alignment_line($aligned_reference);
	}
	foreach my $key (@sorted_keys)
	{
		$output .= _text_alignment_line($aligned_reads->{$key}) . br;
	}	
	foreach my $aligned_reference (@aligned_references)
	{
		$output .= _text_alignment_line($aligned_reference);
	}
	return $output;	
}

sub _strand_char
{
	my ($s) = @_;
	return '<' if ($s == -1);
	return '>' if ($s == +1);
	return '.';
}

sub _text_alignment_line
{
	my ($a, $coords) = @_;
	return $a->{aligned_bases} . "  " . _strand_char($a->{strand}) . "  " . $a->{seq_id} . "\n";
}

sub html_alignment
{
	my $verbose = 0;
	my ($self, $bam_path, $fasta_path, $region, $options) = @_;
	
	my $alignment_info = $self->create_alignment($bam_path, $fasta_path, $region, $options);


	my $output = '';	
	return p . "No reads uniquely align to region." if (!defined $alignment_info);
	$output .= p . "$alignment_info->{message}" if ($alignment_info->{message});
	return $output if (!defined $alignment_info->{aligned_reads});

	my $aligned_reads = $alignment_info->{aligned_reads};
	my @aligned_references = @{$alignment_info->{aligned_references}};
	my $aligned_annotation = $alignment_info->{aligned_annotation};
	my $quality_range = $self->set_quality_range($aligned_reads, $options);
	
	my @sorted_keys = sort { -($aligned_reads->{$a}->{aligned_bases} cmp $aligned_reads->{$b}->{aligned_bases}) } keys %$aligned_reads;

	$output .= style($self->{header_style_string});
	$output .= start_table({-style=>"background-color: rgb(255,255,255)"}) . start_Tr() . start_td({-style=>"font-size:10pt"});
		
	foreach my $aligned_reference (@aligned_references)
	{
		$output .= $self->_html_alignment_line($aligned_reference, 1) . br;
	}
	$output .= $self->_html_alignment_line($aligned_annotation, 0) . br;
		
	foreach my $key (@sorted_keys)
	{
		$output .= $self->_html_alignment_line($aligned_reads->{$key}, 0, $quality_range) . br;
	}	
	$output .= $self->_html_alignment_line($aligned_annotation, 0) . br;
	foreach my $aligned_reference (@aligned_references)
	{
		$output .= $self->_html_alignment_line($aligned_reference, 1) . br;
	}
	$output .= br;

	## create legend information
	
	$output .= start_code . "Base quality scores:&nbsp;" . end_code;
	$output .= $self->_html_alignment_line({aligned_bases => 'ATCG', => aligned_quals => pack('CCCC',0,0,0,0)}, 0,  $quality_range);
	for (my $i=((defined $options->{quality_score_cutoff}) ? 1 : 2); $i<scalar @{$quality_range->{qual_cutoffs}}; $i++)
	{
		my $c = $quality_range->{qual_cutoffs}->[$i];
		$output .= start_code . "&nbsp;&lt;&nbsp;$c&nbsp;&le;&nbsp;" . end_code;
		$output .= $self->_html_alignment_line({aligned_bases => 'ATCG', => aligned_quals => pack('CCCC',$c,$c,$c,$c)}, 0,  $quality_range);
	}
	
	$output .= end_table() . end_Tr() . end_td();	
		
	return $output;
}

sub _html_strand_char
{
	my ($s) = @_;
	return '&lt;' if ($s == -1);
	return '&gt;' if ($s == +1);
	return '.';
}

sub _html_alignment_line
{
	my ($self, $a, $coords, $quality_range) = @_;
	my $output;
	$output .= start_code;
	
	my @split_aligned_bases = split //, $a->{aligned_bases};
	my @split_aligned_quals;
	@split_aligned_quals = split //, $a->{aligned_quals} if ($a->{aligned_quals});
	if (@split_aligned_quals)
	{
		
		if (scalar @split_aligned_bases != scalar @split_aligned_quals)
		
		{
			print "@split_aligned_bases\n";
			print "@split_aligned_quals\n";
		#	$self->throw("unequal aligned base and aligned quals");
		
		}
	}
	
	my $last_color = '';
	for (my $i=0; $i<scalar @split_aligned_bases; $i++)
	{
		my $q = 255; 
		$q = ord($split_aligned_quals[$i]) if (@split_aligned_quals);
		my $b = $split_aligned_bases[$i];	
				
		my $color;		
		## no base quality provided -- assume BEST if not space
		if ($q != 254)
		{
			if (($q == 255) ||  (!defined $quality_range))
			{
				$color = ($b eq ' ') ? 'NC' : "\U$b" . $self->{no_color_index};
			}
			##Note: no color for $q == 254

			elsif ((defined $quality_range) && (!($b =~ m/[.-]/))) #($b =~ m/[NATCGnatcg]/))
			{
				my $color_num = $quality_range->{qual_to_color_index}->[$q];
				#$color = $base_colors_hash->{"\U$b"}->[$color_num];
				$color = "\U$b" . $color_num;
			}
		}
		
		$b = '&nbsp;' if ($b eq ' ');	
		if (not defined $color)
		{
			$color = "UN";
			#older version
			#$color = "color: rgb(0,0,0); background-color: rgb(255,255,255)";
		}

		if ($color ne $last_color)
		{
			$output .= "</font>" if ($last_color);
			$output .= "<font class=\"$color\">";
			$last_color = $color;
		}
		$output .= $b;
		
	}
	$output .= "</font>" if ($last_color);
	$output .= "&nbsp;&nbsp;" . _html_strand_char($a->{strand})  if (defined $a->{strand});
	
	
	##write the seq_id and coords in non-breaking html
	if (defined $a->{seq_id})
	{
		my $seq_id = $a->{seq_id};
		$seq_id =~ s/-/&#8209;/g;
		$seq_id .= "/$a->{start}&#8209;$a->{end}"	if (defined $coords);			
		$output .= "&nbsp;&nbsp;" . $seq_id; 
	}
	
	$output .=  end_code;
	
	return $output;
}


## Workaround to avoid too many open files
my %open_bam_files;

sub create_alignment
{
	my ($self, $bam_path, $fasta_path, $region, $options) = @_;
	my $verbose = $options->{'verbose'};
		
	## Start -- Workaround to avoid too many open files
	if (!defined $open_bam_files{$bam_path.$fasta_path})
	{
		$open_bam_files{$bam_path.$fasta_path} = Bio::DB::Sam->new(-bam =>$bam_path, -fasta=>$fasta_path);
	}
	my $bam = $open_bam_files{$bam_path.$fasta_path};
	## End -- Workaround to avoid too many open files

	#$verbose = 1 if ($region eq "NC_001416‑1:4566-4566");
	my ($seq_id, $start, $end);
	my ($insert_start, $insert_end) = (0, 0);
	#syntax that includes insert counts
	# e.g. NC_001416:4566.1-4566.1
	if ($region =~ m/(.+)\:(\d+)\.(\d+)-(\d+)\.(\d+)/)
	{	
		($seq_id, $start, $insert_start, $end, $insert_end) = ($1, $2, $3, $4, $5);
	}
	elsif ($region =~ m/(.+)\:(\d+)(\.\.|\-)(\d+)/)
	{
		($seq_id, $start, $end) = ($1, $2, $4);
	}
	else
	{
		($seq_id, $start, $end) = split /:|\.\.|\-/, $region;
	}
	my $reference_length = $bam->length($seq_id);
	
	##check the start and end for sanity....	
	$start = 1 if ($start < 1); 
	$end = $reference_length if ($end > $reference_length); 
	return if ($start > $end);
	$region = "$seq_id:$start-$end";
	print "$bam_path  $fasta_path  $region\n" if ($verbose);


	my $aligned_reads;
	my $aligned_annotation;	
	
	my @aligned_references = ( {} );
	## More complex information about multiple reference sequence lines was provided
	if (defined $options->{alignment_reference_info_list})
	{
		@aligned_references = @{$options->{alignment_reference_info_list}};
	}
	foreach my $aligned_reference (@aligned_references)
	{
		$aligned_reference->{seq_id} = $seq_id;
		$aligned_reference->{strand} = 0;
	}
	
	## Need to ignore positions on the left and right that are only overlapped by
	## Redundant reads. Currently Pileup will start and end on those coords.
	## Since it doesn't know about redundancy marking.
	
	my $unique_start;
	my $unique_end;
	
	my $total_reads = 0;
	my $processed_reads = 0;

	## Retrieve all unique alignments overlapping position with "fetch"
	## This lets us know how many slots we need to reserve for alignments.
	my $fetch_function = sub {
		my ($a) = @_;
		#print $a->display_name,' ',$a->cigar_str,"\n";
		my $redundancy = $a->aux_get('X1');		
				
		if ($redundancy == 1)
		{
			$total_reads++;
			return if ($total_reads > $self->{maximum_to_make_alignment});
						
			my $aligned_read;
			$aligned_read->{seq_id} = $a->display_name;
			$aligned_read->{length} = $a->l_qseq;
			$aligned_read->{read_sequence} = $a->qseq;
			$aligned_read->{qual_sequence} = $a->_qscore;

# this biases toward late alignments by how it re-assigns!!			
#			##clear a spot if we're above the limit
#			my $odd_this_one = 1;
#			if ($total_reads > $self->{maximum_to_align})
#			{
#				my @read_keys = keys %$aligned_reads;
#				my $chosen = int(rand($self->{maximum_to_align}+1));
#				
#				##must be a chance of deleting the current one
#				if ($chosen != $self->{maximum_to_align})
#				{
#					delete $aligned_reads->{$read_keys[$chosen]};
#					$aligned_reads->{$a->display_name} = $aligned_read;
#				}
#				else
#				{
#					$odd_this_one = 0;
#				}
#			}

			## save in the hash, creating a spot for each read we will be aligning
			$aligned_reads->{$a->display_name} = $aligned_read;

			## keep track of the earliest and latest coords we see in UNIQUE alignments
			$unique_start = $a->start if (!defined $unique_start || $unique_start > $a->start);
			$unique_end   = $a->end   if (!defined $unique_end   || $unique_end   < $a->end  );				
		}
	};
	$bam->fetch($region, $fetch_function);	
	
	### If there are WAY too many reads, such that a pileup would take forever,
	### then bail.
	if ($total_reads > $self->{maximum_to_make_alignment})
	{
		return { 
			message => "Reads exceeded maximum to display alignment. $total_reads reads. (Limit = $self->{maximum_to_make_alignment})",
		};
	}
	
	my $message;
	if ($total_reads > $self->{maximum_to_align})
	{
		$message = "Only $self->{maximum_to_align} of $total_reads total aligned reads displayed.";
		my $new_aligned_reads;
		my @new_keys = shuffle(keys %$aligned_reads);
		foreach (my $i=0; $i<$self->{maximum_to_align}; $i++)
		{
			$new_aligned_reads->{$new_keys[$i]} = $aligned_reads->{$new_keys[$i]};
		}
		
		$aligned_reads = $new_aligned_reads;
	}
	
	#create the alignment via "pileup"
	my $last_pos;
	my $pileup_function = sub {
    	my ($seq_id,$pos,$pileup) = @_;
		print "POSITION: $pos\n" if ($verbose);

		return if ($pos < $unique_start);
		return if ($pos > $unique_end);

		foreach my $aligned_reference (@aligned_references)
		{
			$aligned_reference->{start} = $pos if (!defined $aligned_reference->{start});
			$aligned_reference->{end} = $pos;
		}		
		
		## Cull the list to those we want to align
		## Other alignments may overlap this position that DO NOT
		## overlap the positions of interest. This removes them.
		@$pileup = grep { defined $aligned_reads->{$_->alignment->display_name} } @$pileup;		

		## Find the maximum indel count
		## We will add gaps to reads with an indel count lower than this.
		my $max_indel = 0;
		my $alignment_spans_position;
		ALIGNMENT: for my $p (@$pileup) 
		{
			## alignment spans the position unless this is the first base...
			$alignment_spans_position->{$p->alignment->display_name} = 1 if ($p->qpos > 0);
			$max_indel = $p->indel if ($p->indel > $max_indel);
		}
		print "MAX INDEL: $max_indel\n" if ($verbose);
		
		## Reference only positions, with no aligned reads
		## are never called, so we keep track of the last position to add them
		if (defined $last_pos && ($last_pos < $pos - 1))
		{
			$last_pos++;
			while ($last_pos < $pos)
			{	
				## READS: add gaps to all
				foreach my $key (keys %$aligned_reads)
				{
					my $aligned_read = $aligned_reads->{$key};
					$aligned_read->{aligned_bases} .= ($alignment_spans_position->{$key}) ? '.' : ' ';
					$aligned_read->{aligned_quals} .= chr(255);					
				}			
			
				## REFERENCE SEQUENCES: add actual bases
				my $ref_bases = $bam->segment($seq_id,$last_pos,$last_pos)->dna;
				foreach my $aligned_reference (@aligned_references)
				{
					my $my_ref_bases = $ref_bases;
					$my_ref_bases = '.' if ($aligned_reference->{truncate_start} && ($last_pos < $aligned_reference->{truncate_start}))
					                 || ($aligned_reference->{truncate_end} && ($last_pos > $aligned_reference->{truncate_end}));
					$aligned_reference->{aligned_bases} .= $my_ref_bases;
					$aligned_reference->{aligned_quals} .= chr(255);
				}
				
				## ANNOTATIONS: add gaps				
				$aligned_annotation->{aligned_bases} .= 
					( (($insert_start == 0) && ($last_pos == $start)) || (($last_pos > $start) && ($last_pos <= $end)) ) 
					? '|' : ' ';
				
				$last_pos++;
			}
		}
		$last_pos = $pos;
		## END adding reference only positions.
				
		## Now add this position to the alignments
		my $updated;
		ALIGNMENT: for my $p (@$pileup) 
		{
			my $a = $p->alignment;
			$updated->{$a->display_name} = 1;
			
			##this setup gives expected behavior for indels!
			my $indel = $p->indel;
			$indel = 0 if ($indel < 0);
			$indel = -1 if ($p->is_del);
			
			## Which read are we on?
			my $aligned_read = $aligned_reads->{$a->display_name};
			$aligned_read->{strand} = ($a->reversed) ? -1 : +1;
			$aligned_read->{reference_start} = $pos if (!defined $aligned_read->{reference_start});
			$aligned_read->{reference_end} = $pos;

			$aligned_read->{start} = $p->qpos+1 if (!defined $aligned_read->{start});
			$aligned_read->{end} = $p->qpos+1;			
			
			## READS: add aligned positions			
			for (my $i=0; $i<=$max_indel; $i++)
			{
				if ($i > $indel)
				{
					print $p->indel . "\n" if ($verbose);
					$aligned_read->{aligned_bases} .= '.';
					$aligned_read->{aligned_quals} .= chr(255);					
				}
				else
				{
					my $quality = $a->qscore->[$p->qpos+$i];
					my $base  = substr($a->qseq, $p->qpos+$i,1);

					if (!$options->{text})
					{
						my $trim_left = $a->aux_get('XL');
						$base = "\L$base" if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
						## alternate coloring scheme
						#$base = chr(ord($base)+128) if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
					  
						my $trim_right = $a->aux_get('XR');	
						$base = "\L$base" if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));
						## alternate coloring scheme					
						#$base = chr(ord($base)+128) if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));						
					}
				
					$aligned_read->{aligned_bases} .= $base;
					$aligned_read->{aligned_quals} .= chr($quality);
				}
				
				print $aligned_read->{aligned_bases} . " " .$aligned_read->{seq_id} . " " . $p->indel . "\n" if ($verbose);
			}
		}
		
		print Dumper($updated) if ($verbose);
		
		## READS: handle those with no 		
		foreach my $key (keys %$aligned_reads)
		{
			next if ($updated->{$key});
			my $aligned_read = $aligned_reads->{$key};
			$aligned_read->{aligned_bases} .= ' ' x ($max_indel+1);
			$aligned_read->{aligned_quals} .= chr(254) x ($max_indel+1);
			
			print $aligned_read->{aligned_bases} . " NOALIGN " . $key . " " . "\n" if ($verbose);
		}
		
		##now handle the reference sequence
		my $ref_base = $bam->segment($seq_id,$pos,$pos)->dna;
		foreach my $aligned_reference (@aligned_references)
		{
			my $my_ref_base = $ref_base;
			$my_ref_base = '.' if ($aligned_reference->{truncate_start} && ($last_pos < $aligned_reference->{truncate_start}))
			                 || ($aligned_reference->{truncate_end} && ($last_pos > $aligned_reference->{truncate_end}));
			
			$aligned_reference->{aligned_bases} .= $my_ref_base;
			$aligned_reference->{aligned_bases} .= '.' x ($max_indel) if ($max_indel > 0);
			$aligned_reference->{aligned_quals} .= chr(255) x ($max_indel+1);
		}
		
		##also update any positions of interest for gaps
		for (my $insert_count=0; $insert_count<=$max_indel; $insert_count++)
		{	
			if (   (($insert_start <= $insert_count) && ($insert_count <= $insert_end) && ($pos == $start) && ($pos == $end))
				|| (($insert_start <= $insert_count) && ($pos == $start) && ($pos != $end)) 
				|| (($pos < $end) && ($pos > $start)) 
				|| (($insert_end <= $insert_count) && ($pos == $end) && ($pos != $start)) )
			{ 
				$aligned_annotation->{aligned_bases} .= '|'; 
			}	
			else
			{
				$aligned_annotation->{aligned_bases} .= ' ';
			}	
		}
	};
	
	$bam->pileup($region, $pileup_function);

    #### IF NOTHING ALIGNED, RETURN undef 
	return undef if (scalar keys %$aligned_reads == 0);

	##now add the unaligned portions of each
	my $max_extend_left = 0;
	my $max_extend_right = 0;
	foreach my $key (keys %$aligned_reads)
	{
		my $aligned_read = $aligned_reads->{$key};
		
		$aligned_read->{aligned_bases} =~ m/^(\s*)\S+(\s*)$/;
		
#		print "\"$aligned_read->{aligned_bases}\"\n";
		
		print "$aligned_read->{start} $aligned_read->{end}\n" if ($verbose);
		
		my $extend_left =  ($aligned_read->{start}-1) - length($1);
		my $extend_right = ($aligned_read->{length}-$aligned_read->{end}) - length($2);

		$max_extend_left = $extend_left if ($extend_left > $max_extend_left);
		$max_extend_right = $extend_right if ($extend_right > $max_extend_right);
	}
	
	print "Extend: $max_extend_left $max_extend_right\n" if ($verbose);
		
	#now add this much to every one
	foreach my $key (keys %$aligned_reads)
	{		
		my $aligned_read = $aligned_reads->{$key};		
		$aligned_read->{aligned_bases} = (' ' x $max_extend_left) . $aligned_read->{aligned_bases} . (' ' x $max_extend_right);
		$aligned_read->{aligned_quals} = (chr(255) x $max_extend_left) . $aligned_read->{aligned_quals} . (chr(255) x $max_extend_right);
	}
			
	###extend reference sequence as requested, be aware of ends of sequence
	## handle left side extending past end
	foreach my $aligned_reference (@aligned_references)
	{	
		my $ref_extend_left = $max_extend_left;
		my $ref_add_left = '';	
		if ($aligned_reference->{start}-$max_extend_left < 1)
		{
			$ref_extend_left = $aligned_reference->{start} - 1;
			$ref_add_left .= '-' x (1 - ($aligned_reference->{start}-$max_extend_left));
		}
		## handle right side extending past end
		my $ref_extend_right = $max_extend_right;
		my $ref_add_right = '';
		if ($aligned_reference->{end}+$max_extend_right > $reference_length)
		{
			$ref_extend_right = $reference_length - $aligned_reference->{end};
			$ref_add_right .= '-' x (($aligned_reference->{end}+$max_extend_right) - $reference_length);
		}
		
		my $pos;
		$pos = $aligned_reference->{start}-$ref_extend_left;
		while ($pos < $aligned_reference->{start})
		{
			my $base = $bam->segment($aligned_reference->{seq_id},$pos,$pos)->dna;	
			$base = '.' if ($aligned_reference->{truncate_start} && ($pos < $aligned_reference->{truncate_start}))
			                 || ($aligned_reference->{truncate_end} && ($pos > $aligned_reference->{truncate_end}));
			$ref_add_left .= $base;
			$pos++;
		}
		
		$pos = $aligned_reference->{end}+1;
		while ($pos <= $aligned_reference->{end}+$ref_extend_right)
		{
			my $base = $bam->segment($aligned_reference->{seq_id},$pos,$pos)->dna;	
			$base = '.' if ($aligned_reference->{truncate_start} && ($pos < $aligned_reference->{truncate_start}))
			                 || ($aligned_reference->{truncate_end} && ($pos > $aligned_reference->{truncate_end}));
			$ref_add_right .= $base;
			$pos++;
		}
		
	#	$ref_add_left .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{start}-$ref_extend_left,$aligned_reference->{start}-1)->dna;	
	#	$ref_add_right .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{end}+1,$aligned_reference->{end}+$ref_extend_right)->dna;
			
		$aligned_reference->{start} -= $ref_extend_left;
		$aligned_reference->{end} += $ref_extend_right;
		$aligned_reference->{aligned_bases} = $ref_add_left . $aligned_reference->{aligned_bases} . $ref_add_right;
		$aligned_reference->{aligned_quals} = (chr(255) x $max_extend_left) . $aligned_reference->{aligned_quals} . (chr(255) x $max_extend_right);
	}
	
	#extend annotation line
	$aligned_annotation->{aligned_bases} = (' ' x $max_extend_left) . $aligned_annotation->{aligned_bases} . (' ' x $max_extend_right);
	
	#now go in and replace the empty space adjacent to each with the rest of the read sequence
	foreach my $key (keys %$aligned_reads)
	{
		my $aligned_read = $aligned_reads->{$key};		
		
		$aligned_read->{aligned_bases} =~ m/^(\s*)(\S+)/;
		my $left_pos =  length($1);
		my $right_pos = $left_pos + length($2);
		
		if ($aligned_read->{start} > 1)
		{
			my $add_seq = substr $aligned_read->{read_sequence}, 0, $aligned_read->{start}-1;
			substr($aligned_read->{aligned_bases}, $left_pos-length($add_seq), length($add_seq)) = "\L$add_seq";

#to color according to quals
#			$add_seq = substr $aligned_read->{qual_sequence}, 0, $aligned_read->{start}-1;

			$add_seq = chr(254) x ($aligned_read->{start}-1);
			substr($aligned_read->{aligned_quals}, $left_pos-length($add_seq), length($add_seq)) = $add_seq;
		}
		if ($aligned_read->{end} < $aligned_read->{length})
		{
			my $add_seq = substr $aligned_read->{read_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
			substr($aligned_read->{aligned_bases}, $right_pos, length($add_seq)) = "\L$add_seq";

#to color according to quals			
#			$add_seq = substr $aligned_read->{qual_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
			
			$add_seq = chr(254) x ($aligned_read->{length}-$aligned_read->{end});
			substr($aligned_read->{aligned_quals}, $right_pos, length($add_seq)) = $add_seq;			
		}
	}
	
	## swap out the ghost seq_ids
	foreach my $ar (@aligned_references)
	{	
		$ar->{seq_id} = $ar->{ghost_seq_id} if (defined $ar->{ghost_seq_id});

		if (defined $ar->{truncate_start})
		{
			$ar->{start} = $ar->{truncate_start};
			my $len = $ar->{end} - $ar->{start} + 1;
			$ar->{start} = $ar->{ghost_start};
			$ar->{end} = $ar->{start} + $ar->{ghost_strand} * ($len - 1);
		} 
		if (defined $ar->{truncate_end})
		{
			$ar->{end} = $ar->{truncate_end};
			my $len = $ar->{end} - $ar->{start} + 1;
			$ar->{end} = $ar->{ghost_end};
			$ar->{start} = $ar->{end} + $ar->{ghost_strand} * ($len - 1);
		}	
	}
	
	### Need to reverse the coords for some
	foreach my $key (keys %$aligned_reads)
	{		
		my $aligned_read = $aligned_reads->{$key};		
		if ($aligned_read->{strand} == -1)
		{
			($aligned_read->{start}, $aligned_read->{end}) = ($aligned_read->{length} - $aligned_read->{start} + 1, $aligned_read->{length} - $aligned_read->{end} + 1);
		}
	}			
			
	return { 
		aligned_reads => $aligned_reads, 
		aligned_references => \@aligned_references, 
		aligned_annotation => $aligned_annotation, 
		message => $message, 
	};
}


=head2 set_quality_range

 Title   : set_quality_range
 Usage   : set_quality_range( min, max );
 Function: 
 Returns : set minimum and maximum qualities
 
=cut

##so these parameters can be set from outside
our $max_quality;
our $min_quality;
our $quality_type;
our @qual_cutoffs;

sub set_quality_range
{
	my ($self, $aligned_reads, $options) = @_;
	
	my @qc;
	my $total = 0;
	my $quality_score_cutoff = 0;
	$quality_score_cutoff = $options->{quality_score_cutoff} if (defined $options->{quality_score_cutoff});	
	
	foreach my $key (keys %$aligned_reads)
	{
		my $aligned_read = $aligned_reads->{$key};
		foreach my $c (split (//, $aligned_read->{qual_sequence}))
		{
			if (ord($c) >= $quality_score_cutoff)
			{
				$qc[ord($c)]++;
				$total++;
			}
		}
	}
	
	my @qual_to_color;
	my @cutoff_percentiles = (0, 0.03, 0.1, 0.3, 0.9, 1.0);
	my $current_cutoff_level = 0;
	
	##set up to this score to the zero level (which is a completely different color)
	my $i;
	for ($i=0; $i<$quality_score_cutoff; $i++)
	{
		$qual_to_color[$i] = $current_cutoff_level;
	}
	$current_cutoff_level++;
	
	my $cumq = 0;	
	while ($i < scalar @qc)
	{
		$cumq += $qc[$i] / $total if (defined $qc[$i]);
		
		#this can increment by at most one per quality score
		if ($cumq > $cutoff_percentiles[$current_cutoff_level])
		{
			$current_cutoff_level++;
		}
		$qual_to_color[$i] = $current_cutoff_level;
		
	} continue {
		$i++;
	}			
		
	#last must be set to max
	$qual_to_color[$i-1] = scalar(@cutoff_percentiles)-1;
	#first must be set to min
	$qual_to_color[$quality_score_cutoff] = 1;	
	
	#if there are at least as many quality scores in existence as
	#there are color levels to assign....
	if ((scalar(@qual_to_color) > scalar(@cutoff_percentiles)-1))
	{
		#...redistribute such that there are no jumps in quality level
		my $gap = 1;
		while ($gap)
		{
			$gap = 0;
			my $last = 0;
			for (my $i=0; $i<scalar @qual_to_color; $i++)
			{
				if ($qual_to_color[$i] > $last + 1)
				{
					$qual_to_color[$i-1]++;
					$gap = 1;
				}
				$last = $qual_to_color[$i];
			}
		}
	}			
		
	##finally, this sets the cutoff levels
	my $last = 0;
	my @cutoff_levels = ($quality_score_cutoff);
	for (my $i=$quality_score_cutoff; $i<scalar @qual_to_color; $i++)
	{
		if ($qual_to_color[$i] > $last)
		{
			push @cutoff_levels, $i;
			$last = $qual_to_color[$i];
		}
	}	

	return {qual_to_color_index => \@qual_to_color, qual_cutoffs => \@cutoff_levels };
}

return 1;

