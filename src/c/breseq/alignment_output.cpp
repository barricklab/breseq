/*****************************************************************************

AUTHORS

	Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
	David B. Knoester

LICENSE AND COPYRIGHT

	Copyright (c) 2008-2010 Michigan State University
	Copyright (c) 2011 The University of Texas at Austin

	breseq is free software; you can redistribute it and/or modify it under the
	terms the GNU General Public License as published by the Free Software
	Foundation; either version 1, or (at your option) any later version.

 *****************************************************************************/

#include "breseq/alignment_output.h"
#include "backup/alignment_output.h"
#include "breseq/common.h"
#include "breseq/pileup.h"
#include "assert.h"
#include "boost/optional.hpp"

using namespace std;

namespace breseq {

	bool verbose = true; //TODO Boost::options
	bool debug = true;

	alignment_output_pileup::alignment_output_pileup(const string& bam, const string& fasta, const uint32_t maximum_to_align)
	: pileup_base(bam, fasta)
	, maximum_to_align(maximum_to_align) {
	}

	alignment_output_pileup::~alignment_output_pileup() {
	}

	alignment_output::alignment_output(string bam, string fasta, uint32_t maximum_to_align)
	: m_alignment_output_pileup_object(bam, fasta, maximum_to_align) {
	}

	void alignment_output::create_alignment(const string bam, const string fasta, const string region) {
		//sub create_alignment
		//{
		//	my ($self, $bam_path, $fasta_path, $region, $options) = @_;
		//	my $verbose = $options->{'verbose'};
		//
		//	my ($seq_id, $start, $end, $insert_start, $insert_end) = Breseq::Shared::region_to_coords($region, $bam);
		//	$region = "$seq_id:$start-$end";
		//	print "$bam_path  $fasta_path  $region\n" if ($verbose);
		//
		//	my $reference_length = $bam->length($seq_id);
		//	my $aligned_reads;
		//	my $aligned_annotation;
		//
		//	my @aligned_references = ( {} );
		//	////## More complex information about multiple reference sequence lines was provided
		//	////if (defined $options->{alignment_reference_info_list})
		//	////{
		//	////	@aligned_references = @{$options->{alignment_reference_info_list}};
		//	////}
		//	////foreach my $aligned_reference (@aligned_references)
		//	////{
		//	////	$aligned_reference->{seq_id} = $seq_id;
		//	////	$aligned_reference->{strand} = 0;
		//	////}
		//
		//	## Need to ignore positions on the left and right that are only overlapped by
		//	## Redundant reads. Currently Pileup will start and end on those coords.
		//	## Since it doesn't know about redundancy marking.
		//
		//	my $unique_start;
		//	my $unique_end;
		//
		//	my $total_reads = 0;
		//	my $processed_reads = 0;
		//
		m_alignment_output_pileup_object.unique_start = 0; //TODO Move to constructor
		m_alignment_output_pileup_object.unique_end = 0; //TODO Move to constructor
		m_alignment_output_pileup_object.total_reads = 0; //TODO Move to constructor
		m_alignment_output_pileup_object.processed_reads =0; //TODO Move to constructor

		//*Call
		//* do fetch_callback
		//	$bam->fetch($region, $fetch_function);

		//
		//	### There are not uniquely aligned reads...
		//	if (!defined $unique_start || !defined $unique_end)
		//	{
		//		return;
		//	}

		//
		//	### If there are WAY too many reads, such that a pileup might take forever, bail...
		//	if ($self->{maximum_to_make_alignment} && ($total_reads > $self->{maximum_to_make_alignment}))
		//	{
		//		return {
		//			message => "Reads exceeded maximum to display alignment. $total_reads reads. (Limit = $self->{maximum_to_make_alignment})",
		//		};
		//	}
		//
			m_alignment_output_pileup_object.do_fetch(region);

			if(debug)
			{
			  cout << "End do_fetch function, variables =";
				cout << "unique_start: " << m_alignment_output_pileup_object.unique_start << "\t";
				cout << "unique_end: " << m_alignment_output_pileup_object.unique_end << endl;
			}

			if ((m_alignment_output_pileup_object.unique_start == 0) || (m_alignment_output_pileup_object.unique_end == 0))
			{
			  cout << "No unique start or end initialized" << endl;
				return;
			}

      if (m_alignment_output_pileup_object.total_reads > m_alignment_output_pileup_object.maximum_to_align)
			{
				cout << "Reads exceeded maximum to display alignment. ";
				cout << m_alignment_output_pileup_object.total_reads << " reads. ";
				cout << "(Limit = " << m_alignment_output_pileup_object.maximum_to_align;
				cout << ")" << endl;
				return;
			}
		//	my $message;
		//////	if ($self->{maximum_to_align} && ($total_reads > $self->{maximum_to_align}))
		//////	{
		//////		$message = "Only $self->{maximum_to_align} of $total_reads total aligned reads displayed.";
		//////		my $new_aligned_reads;
		//////		my @new_keys = shuffle(keys %$aligned_reads);
		//////		foreach (my $i=0; $i<$self->{maximum_to_align}; $i++)
		//////		{
		//////			$new_aligned_reads->{$new_keys[$i]} = $aligned_reads->{$new_keys[$i]};
		//////		}
		//
		//////		$aligned_reads = $new_aligned_reads;
		//////	}
		//

		// * Call
		// * do_pileup ();
		// *
		//
		//	$bam->pileup($region, $pileup_function);
 	m_alignment_output_pileup_object.do_pileup(region);


		//
		//    #### IF NOTHING ALIGNED, RETURN undef
		//	return undef if (scalar keys %$aligned_reads == 0);
		//
		//	##now add the unaligned portions of each
		//	my $max_extend_left = 0;
		//	my $max_extend_right = 0;
		//	foreach my $key (keys %$aligned_reads)
		//	{
		//		my $aligned_read = $aligned_reads->{$key};
		//
		//		$aligned_read->{aligned_bases} =~ m/^(\s*)\S+(\s*)$/;
		//
		//#		print "\"$aligned_read->{aligned_bases}\"\n";
		//
		//		print "$aligned_read->{start} $aligned_read->{end}\n" if ($verbose);
		//
		//		my $extend_left =  ($aligned_read->{start}-1) - length($1);
		//		my $extend_right = ($aligned_read->{length}-$aligned_read->{end}) - length($2);
		//
		//		$max_extend_left = $extend_left if ($extend_left > $max_extend_left);
		//		$max_extend_right = $extend_right if ($extend_right > $max_extend_right);
		//	}
		//
		//	print "Extend: $max_extend_left $max_extend_right\n" if ($verbose);
		//
		//	#now add this much to every one
		//	foreach my $key (keys %$aligned_reads)
		//	{
		//		my $aligned_read = $aligned_reads->{$key};
		//		$aligned_read->{aligned_bases} = (' ' x $max_extend_left) . $aligned_read->{aligned_bases} . (' ' x $max_extend_right);
		//		$aligned_read->{aligned_quals} = (chr(255) x $max_extend_left) . $aligned_read->{aligned_quals} . (chr(255) x $max_extend_right);
		//	}
		//
		//	###extend reference sequence as requested, be aware of ends of sequence
		//	## handle left side extending past end
		//	foreach my $aligned_reference (@aligned_references)
		//	{
		//		my $ref_extend_left = $max_extend_left;
		//		my $ref_add_left = '';
		//		if ($aligned_reference->{start}-$max_extend_left < 1)
		//		{
		//			$ref_extend_left = $aligned_reference->{start} - 1;
		//			$ref_add_left .= '-' x (1 - ($aligned_reference->{start}-$max_extend_left));
		//		}
		//		## handle right side extending past end
		//		my $ref_extend_right = $max_extend_right;
		//		my $ref_add_right = '';
		//		if ($aligned_reference->{end}+$max_extend_right > $reference_length)
		//		{
		//			$ref_extend_right = $reference_length - $aligned_reference->{end};
		//			$ref_add_right .= '-' x (($aligned_reference->{end}+$max_extend_right) - $reference_length);
		//		}
		//
		//		my $pos;
		//		$pos = $aligned_reference->{start}-$ref_extend_left;
		//		while ($pos < $aligned_reference->{start})
		//		{
		//			my $base = $bam->segment($aligned_reference->{seq_id},$pos,$pos)->dna;
		//			$base = '.' if ($aligned_reference->{truncate_start} && ($pos < $aligned_reference->{truncate_start}))
		//			                 || ($aligned_reference->{truncate_end} && ($pos > $aligned_reference->{truncate_end}));
		//			$ref_add_left .= $base;
		//			$pos++;
		//		}
		//
		//		$pos = $aligned_reference->{end}+1;
		//		while ($pos <= $aligned_reference->{end}+$ref_extend_right)
		//		{
		//			my $base = $bam->segment($aligned_reference->{seq_id},$pos,$pos)->dna;
		//			$base = '.' if ($aligned_reference->{truncate_start} && ($pos < $aligned_reference->{truncate_start}))
		//			                 || ($aligned_reference->{truncate_end} && ($pos > $aligned_reference->{truncate_end}));
		//			$ref_add_right .= $base;
		//			$pos++;
		//		}
		//
		//	#	$ref_add_left .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{start}-$ref_extend_left,$aligned_reference->{start}-1)->dna;
		//	#	$ref_add_right .= $bam->segment($aligned_reference->{seq_id},$aligned_reference->{end}+1,$aligned_reference->{end}+$ref_extend_right)->dna;
		//
		//		$aligned_reference->{start} -= $ref_extend_left;
		//		$aligned_reference->{end} += $ref_extend_right;
		//		$aligned_reference->{aligned_bases} = $ref_add_left . $aligned_reference->{aligned_bases} . $ref_add_right;
		//		$aligned_reference->{aligned_quals} = (chr(255) x $max_extend_left) . $aligned_reference->{aligned_quals} . (chr(255) x $max_extend_right);
		//	}
		//
		//	#extend annotation line
		//	$aligned_annotation->{aligned_bases} = (' ' x $max_extend_left) . $aligned_annotation->{aligned_bases} . (' ' x $max_extend_right);
		//
		//	#now go in and replace the empty space adjacent to each with the rest of the read sequence
		//	foreach my $key (keys %$aligned_reads)
		//	{
		//		my $aligned_read = $aligned_reads->{$key};
		//
		//		$aligned_read->{aligned_bases} =~ m/^(\s*)(\S+)/;
		// *
		// *  FIND LENGTH OF THE NON-SPACE PART OF READ
		//		my $left_pos =  length($1);
		//		my $right_pos = $left_pos + length($2);
		//
		//		if ($aligned_read->{start} > 1)
		//		{
		//			my $add_seq = substr $aligned_read->{read_sequence}, 0, $aligned_read->{start}-1;
		//			substr($aligned_read->{aligned_bases}, $left_pos-length($add_seq), length($add_seq)) = "\L$add_seq";
		//
		//#to color according to quals
		//#			$add_seq = substr $aligned_read->{qual_sequence}, 0, $aligned_read->{start}-1;
		//
		//			$add_seq = chr(254) x ($aligned_read->{start}-1);
		//			substr($aligned_read->{aligned_quals}, $left_pos-length($add_seq), length($add_seq)) = $add_seq;
		//		}
		//		if ($aligned_read->{end} < $aligned_read->{length})
		//		{
		//			my $add_seq = substr $aligned_read->{read_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
		//			substr($aligned_read->{aligned_bases}, $right_pos, length($add_seq)) = "\L$add_seq";
		//
		//#to color according to quals
		//#			$add_seq = substr $aligned_read->{qual_sequence}, $aligned_read->{end}, $aligned_read->{length}-$aligned_read->{end};
		//
		//			$add_seq = chr(254) x ($aligned_read->{length}-$aligned_read->{end});
		//			substr($aligned_read->{aligned_quals}, $right_pos, length($add_seq)) = $add_seq;
		//		}
		//	}
		//
		//	## swap out the ghost seq_ids
		//////	foreach my $ar (@aligned_references)
		//////	{
		//		$ar->{seq_id} = $ar->{ghost_seq_id} if (defined $ar->{ghost_seq_id});
		//
		//		if (defined $ar->{truncate_start})
		//		{
		//			$ar->{start} = $ar->{truncate_start};
		//			my $len = $ar->{end} - $ar->{start} + 1;
		//			$ar->{start} = $ar->{ghost_start};
		//			$ar->{end} = $ar->{start} + $ar->{ghost_strand} * ($len - 1);
		//		}
		//		if (defined $ar->{truncate_end})
		//		{
		//			$ar->{end} = $ar->{truncate_end};
		//			my $len = $ar->{end} - $ar->{start} + 1;
		//			$ar->{end} = $ar->{ghost_end};
		//			$ar->{start} = $ar->{end} + $ar->{ghost_strand} * ($len - 1);
		//		}
		//////	}
		//
		//	### Need to reverse the coords for some
		//	foreach my $key (keys %$aligned_reads)
		//	{
		//		my $aligned_read = $aligned_reads->{$key};
		//		if ($aligned_read->{strand} == -1)
		//		{
		//			($aligned_read->{start}, $aligned_read->{end}) = ($aligned_read->{length} - $aligned_read->{start} + 1, $aligned_read->{length} - $aligned_read->{end} + 1);
		//		}
		//	}
		//
		//	return {
		//		aligned_reads => $aligned_reads,
		//		aligned_references => \@aligned_references,
		//		aligned_annotation => $aligned_annotation,
		//		message => $message,
		//	};
		//
	} //End create alignment

	string alignment_output::html_alignment(const string region) {

		string s;

	
	//	m_alignment_output_object.do_pileup(region);

		return s;
	}

	/*! Called for each position.
	 */
void alignment_output_pileup::pileup_callback(const pileup& p) {
		//	#create the alignment via "pileup"
		//	my $last_pos;
		//	my $pileup_function = sub {
		//    	my ($seq_id,$pos,$pileup) = @_;
		//		print "POSITION: $pos\n" if ($verbose);
		//
		//		return if ($pos < $unique_start);
		//		return if ($pos > $unique_end);
		//
		//		foreach my $aligned_reference (@aligned_references)
		//		{
		//			$aligned_reference->{start} = $pos if (!defined $aligned_reference->{start});
		//			$aligned_reference->{end} = $pos;
		//		}
		//
		uint32_t last_pos = 0; // TODO move to constructor
		uint32_t pos = p.position_1();
		if (verbose) //TODO Boost::options
			cout << "POSITION: " << pos << endl;
		
	
		if ((pos < unique_start) || (pos > unique_end))
			return;

		for (uint32_t index_aligned_reference;		
					 index_aligned_reference < aligned_references.size();
					 index_aligned_reference++)
		{
			aligned_references[index_aligned_reference].start = pos; //TODO check if start is defined
			aligned_references[index_aligned_reference].end =pos;
		}//?Nothing gets defined for single given region case?

		//		## Cull the list to those we want to align
		//		## Other alignments may overlap this position that DO NOT
		//		## overlap the positions of interest. This removes them.
		//		@$pileup = grep { defined $aligned_reads->{$_->alignment->display_name} } @$pileup;
		//
		//		## Find the maximum indel count
		//		## We will add gaps to reads with an indel count lower than this.
		//		my $max_indel = 0;
		//		my $alignment_spans_position;
    //
		//		ALIGNMENT: for my $p (@$pileup)
		//		{
		//			## alignment spans the position unless this is the first base...
		//			$alignment_spans_position->{$p->alignment->display_name} = 1 if ($p->qpos > 0);
		//			$max_indel = $p->indel if ($p->indel > $max_indel);
		//		}
		//		print "MAX INDEL: $max_indel\n" if ($verbose);
		//
	 max_indel = 0; //TODO move to Constructor
	 for (pileup :: const_iterator itr_pileup = p.begin();
	  		  itr_pileup != p.end() ; itr_pileup ++) {
	   if (itr_pileup->indel() > max_indel)
		 {
		   max_indel = itr_pileup->indel();
		 }
		 alignment_spans_position[itr_pileup->query_name()]
			 =itr_pileup->is_alignment_spanning_position();
	 }
   if (verbose) //TODO Boost::options
	   cout << "MAX INDEL: " << max_indel << endl;


		//		## Reference only positions, with no aligned reads
		//		## are never called, so we keep track of the last position to add them
		//		if (defined $last_pos && ($last_pos < $pos - 1))
		//		{
		//			$last_pos++;
		//			while ($last_pos < $pos)
		//			{
		//				## READS: add gaps to all
		//				foreach my $key (keys %$aligned_reads)
		//				{
		//					my $aligned_read = $aligned_reads->{$key};
		//					$aligned_read->{aligned_bases} .= ($alignment_spans_position->{$key}) ? '.' : ' ';
		//					$aligned_read->{aligned_quals} .= chr(255);
		//				}
		//
		//				## REFERENCE SEQUENCES: add actual bases
		//				my $ref_bases = $bam->segment($seq_id,$last_pos,$last_pos)->dna;
		//				foreach my $aligned_reference (@aligned_references)
		//				{
		//					my $my_ref_bases = $ref_bases;
		//////					$my_ref_bases = '.' if ($aligned_reference->{truncate_start} && ($last_pos < $aligned_reference->{truncate_start}))
		//////					                 || ($aligned_reference->{truncate_end} && ($last_pos > $aligned_reference->{truncate_end}));
		//					$aligned_reference->{aligned_bases} .= $my_ref_bases;
		//					$aligned_reference->{aligned_quals} .= chr(255);
		//				}
		//
		//				## ANNOTATIONS: add gaps
		//				$aligned_annotation->{aligned_bases} .=
		//					( (($insert_start == 0) && ($last_pos == $start)) || (($last_pos > $start) && ($last_pos <= $end)) )
		//					? '|' : ' ';
		//
		//				$last_pos++;
		//			}
		//		}
	 if( ( last_pos != 0 ) && ( last_pos < pos) )
	{
		  last_pos++;
				while (last_pos < p.position_1())
				{
					for (map<string, struct_aligned_read>::const_iterator itr_reads = aligned_reads.begin(); //TODO typedef map
								itr_reads != aligned_reads.end(); itr_reads++) {
					  aligned_reads[itr_reads->first].aligned_read_bases += ".";
						//START //Unique
					}
				 	}
		}
	last_pos = pos;



		//		$last_pos = $pos;
		//		## END adding reference only positions.
		//
		//		## Now add this position to the alignments
		//		my $updated;
		//		ALIGNMENT: for my $p (@$pileup)
		//		{
		//			my $a = $p->alignment;
		//			$updated->{$a->display_name} = 1;
		//
		//			##this setup gives expected behavior for indels!
		//			my $indel = $p->indel;
		//			$indel = 0 if ($indel < 0);
		//			$indel = -1 if ($p->is_del);
		//
		//			## Which read are we on?
		//			my $aligned_read = $aligned_reads->{$a->display_name};
		//			$aligned_read->{strand} = ($a->reversed) ? -1 : +1;
		//			$aligned_read->{reference_start} = $pos if (!defined $aligned_read->{reference_start});
		//			$aligned_read->{reference_end} = $pos;
		//
		//			$aligned_read->{start} = $p->qpos+1 if (!defined $aligned_read->{start});
		//			$aligned_read->{end} = $p->qpos+1;
		//
		//			## READS: add aligned positions
		//			for (my $i=0; $i<=$max_indel; $i++)
		//			{
		//				if ($i > $indel)
		//				{
		//					print $p->indel . "\n" if ($verbose);
		//					$aligned_read->{aligned_bases} .= '.';
		//					$aligned_read->{aligned_quals} .= chr(255);
		//				}
		//				else
		//				{
		//					my $quality = $a->qscore->[$p->qpos+$i];
		//					my $base  = substr($a->qseq, $p->qpos+$i,1);
		//
		//					if (!$options->{text})
		//					{
		//						my $trim_left = $a->aux_get('XL');
		//						$base = "\L$base" if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
		//						## alternate coloring scheme
		//						#$base = chr(ord($base)+128) if ((defined $trim_left) && ($p->qpos+1 <= $trim_left));
		//
		//						my $trim_right = $a->aux_get('XR');
		//						$base = "\L$base" if ((defined $trim_right) && ($a->l_qseq-$p->qpos <= $trim_right));
		//						## alternate coloring scheme
		//						#$base = chr(ord($base)+128) if ((defined $trim_right) && ($a->query->length-$p->qpos <= $trim_right));
		//					}
		//
		//					$aligned_read->{aligned_bases} .= $base;
		//					$aligned_read->{aligned_quals} .= chr($quality);
		//				}
		//
		//				print $aligned_read->{aligned_bases} . " " .$aligned_read->{seq_id} . " " . $p->indel . "\n" if ($verbose);
		//			}
		//		}
		//
		//		print Dumper($updated) if ($verbose);
		//
		//		## READS: handle those with no
		//		foreach my $key (keys %$aligned_reads)
		//		{
		//			next if ($updated->{$key});
		//			my $aligned_read = $aligned_reads->{$key};
		//			$aligned_read->{aligned_bases} .= ' ' x ($max_indel+1);
		//			$aligned_read->{aligned_quals} .= chr(254) x ($max_indel+1);
		//
		//			print $aligned_read->{aligned_bases} . " NOALIGN " . $key . " " . "\n" if ($verbose);
		//		}
		//
		//		##now handle the reference sequence
		//		my $ref_base = $bam->segment($seq_id,$pos,$pos)->dna;
		//		foreach my $aligned_reference (@aligned_references)
		//		{
		//			my $my_ref_base = $ref_base;
		//			$my_ref_base = '.' if ($aligned_reference->{truncate_start} && ($last_pos < $aligned_reference->{truncate_start}))
		//			                 || ($aligned_reference->{truncate_end} && ($last_pos > $aligned_reference->{truncate_end}));
		//
		//			$aligned_reference->{aligned_bases} .= $my_ref_base;
		//			$aligned_reference->{aligned_bases} .= '.' x ($max_indel) if ($max_indel > 0);
		//			$aligned_reference->{aligned_quals} .= chr(255) x ($max_indel+1);
		//		}
		//
		//		##also update any positions of interest for gaps
		//		for (my $insert_count=0; $insert_count<=$max_indel; $insert_count++)
		//		{
		//			if (   (($insert_start <= $insert_count) && ($insert_count <= $insert_end) && ($pos == $start) && ($pos == $end))
		//				|| (($insert_start <= $insert_count) && ($pos == $start) && ($pos != $end))
		//				|| (($pos < $end) && ($pos > $start))
		//				|| (($insert_end <= $insert_count) && ($pos == $end) && ($pos != $start)) )
		//			{
		//				$aligned_annotation->{aligned_bases} .= '|';
		//			}
		//			else
		//			{
		//				$aligned_annotation->{aligned_bases} .= ' ';
		//			}
		//		}
		//

	}

	/*! Called for each read alignment.
	 */
	void alignment_output_pileup::fetch_callback(const alignment& a) {

		//
		//	## Retrieve all unique alignments overlapping position with "fetch"
		//	## This lets us know how many slots we need to reserve for alignments.
		//	my $fetch_function = sub {
		//		my ($a) = @_;
		//		#print $a->display_name,' ',$a->cigar_str,"\n";
		//		my $redundancy = $a->aux_get('X1');
		//
		//		if ((!defined $redundancy) || ($redundancy == 1))
		//		{
		//			$total_reads++;
		//
		//			return undef if ($self->{maximum_to_make_alignment} && ($total_reads > $self->{maximum_to_make_alignment}));
		//
		//			my $aligned_read;
		//			$aligned_read->{seq_id} = $a->display_name;
		//			$aligned_read->{length} = $a->l_qseq;
		//			$aligned_read->{read_sequence} = $a->qseq;
		//			$aligned_read->{qual_sequence} = $a->_qscore;
		//
		//			## save in the hash, creating a spot for each read we will be aligning
		//			$aligned_reads->{$a->display_name} = $aligned_read;
		//
		//			## keep track of the earliest and latest coords we see in UNIQUE alignments
		//			$unique_start = $a->start if (!defined $unique_start || $unique_start > $a->start);
		//			$unique_end   = $a->end   if (!defined $unique_end   || $unique_end   < $a->end  );
		//		}
		//	};
		if (!a.is_redundant())
		{
			total_reads++;

			if(maximum_to_align < total_reads)
				return;

			struct_aligned_read aligned_read;
			aligned_read.seq_id = a.query_name();
			aligned_read.length = a.query_length();
			aligned_read.read_sequence = a.query_char_sequence();
			aligned_read.qual_sequence = string((char*)a.quality_scores());

			aligned_reads[aligned_read.seq_id] = aligned_read;

			if((unique_start == 0) || (unique_start > a.reference_start_1()))
			{
			  unique_start = a.reference_start_1();
			}
			if((unique_end == 0) || (unique_end < a.reference_end_1()))
			{
				unique_end = a.reference_end_1();
			}
			if(debug)
			{
				cout << "unique_start : " << unique_start;
			  cout << "\t" << "unique_end: " << unique_end;
				cout << "\t" << "seq_id: " << aligned_read.seq_id;
				cout << "\t" << "query_start: " << a.query_start_1();
				cout << "\t" << "query_end: " << a.query_end_1();
				cout << endl;
			}
			}
	}


	

} // namespace breseq

