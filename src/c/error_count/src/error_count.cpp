#include <iostream>
#include <sam.h>
#include "error_count.h"


/*! User-defined struct that holds information to be used by the pileup function.
 */
struct user_data {
	//! Constructor.
	user_data() : in(0), unique_only_coverage(0) { }

	samfile_t *in; //!< BAM file handle.
	int unique_only_coverage; //!< ??
};


/*! Pileup callback for error count.
 
 \param tid chromosome id (?)
 \param pos starting position of this alignment
 \param n number of alignments at this position
 \param pl array of n alignments occuring at this position
 */
int error_count_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
	using namespace std;
	user_data* ud = reinterpret_cast<user_data*>(data);
	
	if((pos % 10000) == 0) {
		cerr << "    POSITION:" << pos << endl;
	}
	
	// for each alignment within this pileup:
	for(int i=0; i<n; ++i) {
		const bam_pileup1_t* a=&pl[i];
		
		// is this is a redundant alignment?
		int32_t redundant=bam_aux2i(bam_aux_get(a->b,"X1"));
		if(redundant != 1) {
			continue;	// no; don't process non-unique reads
		}

		// unique_only_position = 0;
		++ud->unique_only_coverage;
		
		//		
		
		//		my $indel = $p->indel;
		//		$indel = 0 if ($indel < 0);
		//		$indel = -1 if ($p->is_del);
		
		//		my $qseq = $a->qseq;
		//		my $qpos = $p->qpos;
		//		my $qscore = $a->qscore;
		//		my $reversed = $a->reversed;
		//		my $query_end = $a->query->end;
		//		my $query_start = $a->query->start;
		//		
		//		my $fastq_file_index = $a->aux_get('X2');
		
			}
	
	/*
	 ## populated by pileup

	 my $unique_only_coverage;
	 my $ref_seq_string = $ref_seq_info->{ref_strings}->{$seq_id};
	 my $com_seq_string = $ref_seq_string;
	 $com_seq_string =~ tr/ATCG/TAGC/;
	 
	 my $pileup_function = 
	 
	 sub {
	 my ($seq_id,$pos,$pileup) = @_;
	 
	 
	 my @ref_base; #index is 'reversed'
	 $ref_base[0] = substr $ref_seq_string, $pos-1, 1;
	 $ref_base[1] = substr $com_seq_string, $pos-1, 1;
	 #$ref_base[0] = $bam->segment($seq_id,$pos,$pos)->dna;
	 #$ref_base[1] = $bam->segment($seq_id,$pos,$pos)->seq->revcom->seq;
	 
	 my $unique_only_position = 1;
	 my $unique_coverage = 0;
	 
	 
	 ALIGNMENT: for my $p (@$pileup) {
	 my $a = $p->alignment;
	 
	 //				##this setup gives expected behavior from indel!
	 my $indel = $p->indel;
	 $indel = 0 if ($indel < 0);
	 $indel = -1 if ($p->is_del);
	 
	 ##
	 # In all that follows, be sure to keep track of strandedness of mutations!
	 ##
	 
	 # (1) base substitutions
	 #     e.g. 'AG' key for observing a G in read at a place where the reference has an A
	 #     IMPROVE by keeping track of flanking base context and scores?
	 #     this would, for example, penalize low scoring sequences more
	 
	 # (2) deletion in read relative to reference
	 #     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
	 #     how does one give a quality score? Use quality score of the next base in
	 #     the read, i.e. where the deleted base would have been in the read
	 #     PROBLEM what do we do about multiple base deletions?
	 
	 # (3) insertion in read relative to reference
	 #     e.g. '.A' key for observing an A in a read at a position where the reference has no base
	 #     how does one give a quality score? - 
	 #     for reference observations: average the quality scores of the surrounding bases in the read (round down)
	 #     for mutation observations: the quality score of the inserted base
	 
	 ## For (2) and (3) train error model only on single base insertions or deletions.
	 
	 ## taking into account neighborhood quality doesn't seem to improve calculations
	 # my $neighborhood_quality = 0;				
	 # my $start_neighborhood = $qpos - 2;
	 # $start_neighborhood = $query_start-1 if ($start_neighborhood < $query_start-1);
	 # my $end_neighborhood = $qpos + 2;
	 # $end_neighborhood = $query_end-1 if ($end_neighborhood > $query_end-1);
	 # foreach my $i ($start_neighborhood..$end_neighborhood) 
	 # {
	 # 	$neighborhood_quality += $qscore->[$i];
	 # }
	 # $neighborhood_quality = POSIX::floor(($neighborhood_quality / ($end_neighborhood - $start_neighborhood + 1)) / 10);
	 
	 //# (1) base substitutions
	 if (!$p->is_del) {
	 my $base  = substr($qseq, $qpos,1);
	 next if ($base =~ /[nN]/);
	 
	 my $quality = $qscore->[$qpos];
	 $base = FastqLite::revcom($base) if ($reversed);
	 my $key = $ref_base[$reversed] . $base; 
	 $error_hash->[$fastq_file_index]->{$quality}->{$key}++;
	 #$complex_error_hash->[$fastq_file_index]->{$neighborhood_quality}->{$quality}->{$key}++;
	 
	 //## also add an observation of a non-gap non-gap
	 if ($qpos+1 < $query_end)
	 {
	 my $next_quality = $qscore->[$qpos+1];
	 my $avg_quality = POSIX::floor( ($quality + $next_quality) / 2);
	 $error_hash->[$fastq_file_index]->{$avg_quality}->{'..'}++;
	 #$complex_error_hash->[$fastq_file_index]->{$neighborhood_quality}->{$avg_quality}->{'..'}++;
	 
	 }
	 }
	 
	 //# (2) deletion in read relative to reference
	 //#     -- only count if there is no indel after this posision
	 elsif ($p->indel == 0)
	 {
	 my $quality = $qscore->[$qpos+(1-$reversed)];
	 #print $a->qname . " " . $pos . " " . $ref_base[$a->reversed] . " " . $quality . "\n";
	 my $key = $ref_base[$reversed] . '.'; 
	 $error_hash->[$fastq_file_index]->{$quality}->{$key}++;	
	 #$complex_error_hash->[$fastq_file_index]->{$neighborhood_quality}->{$quality}->{$key}++;						
	 
	 }
	 
	 //# (3) insertion in read relative to reference
	 //#     -- at the next position in the read
	 //#     -- only count if an indel = +1, meaning a single-base insertion
	 if ($p->indel == +1)
	 {
	 my $base  = substr($qseq,$qpos+1,1);
	 next if ($base =~ /[nN]/);
	 my $quality = $a->qscore->[$qpos+1];
	 #print $a->qname . " " . $pos . " " . $ref_base[$a->reversed] . " " . $quality . "\n";
	 my $key = '.' . $base; 
	 $error_hash->[$fastq_file_index]->{$quality}->{$key}++;
	 #$complex_error_hash->[$fastq_file_index]->{$neighborhood_quality}->{$quality}->{$key}++;
	 }	
	 
	 } #end ALIGNMENT
	 
	 //# record unique only coverage
	 $unique_only_coverage->[$unique_coverage]++ if ($unique_only_position);
	 }; #end $pileup_function
	 */	
	
	return 0;  
}


/*! Calculate the errors in the given BAM file.
 */
void breseq::error_count(const std::string& bam) {
	using namespace std;
	
	user_data ud;
	ud.in = samopen(bam.c_str(), "rb", 0);

	// for each sequence id...
	for(int i=0; i<ud.in->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << ud.in->header->target_name[i] << endl;
		cerr << "  LENGTH: " << ud.in->header->target_len[i] << endl;
		// pileup?
	}
	
	sampileup(ud.in, -1, error_count_callback, &ud);
	samclose(ud.in);
}

/*
 sub count
 {
 my ($settings, $summary, $ref_seq_info) = @_;
 
 my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
 my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
 my $bam = Bio::DB::Sam->new(-fasta => $reference_fasta_file_name, -bam => $reference_bam_file_name);
 my @seq_ids = $bam->seq_ids;
 
 ## populated by pileup
 my $error_hash = [];		#list by fastq file index
 #my $complex_error_hash = [];		#list by fastq file index
 
 foreach my $seq_id (@seq_ids)
 {							
 my $sequence_length = $bam->length($seq_id);
 print STDERR "  REFERENCE: $seq_id\n";
 print STDERR "  LENGTH: $sequence_length\n";
 
 ## populated by pileup
 my $unique_only_coverage;
 my $ref_seq_string = $ref_seq_info->{ref_strings}->{$seq_id};
 my $com_seq_string = $ref_seq_string;
 $com_seq_string =~ tr/ATCG/TAGC/;
 
 // pileup function here.
 
 $bam->pileup($seq_id,$pileup_function);
 #$bam->pileup("$seq_id:1-10000", $pileup_function);
 
 ## save the unique only coverage distribution
 my $this_unique_only_coverage_distribution_file_name = $settings->file_name('unique_only_coverage_distribution_file_name', {'@'=>$seq_id});
 save_unique_coverage_distribution_file($this_unique_only_coverage_distribution_file_name, $unique_only_coverage);
 }
 
 foreach my $read_file ($settings->read_files)
 {
 my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
 my $this_error_counts_file_name = $settings->file_name('error_counts_file_name', {'#'=> $read_file});
 save_error_file($this_error_counts_file_name, $error_hash->[$fastq_file_index]);
 }
 
 # foreach my $read_file ($settings->read_files)
 # {
 # 	my $fastq_file_index = $settings->read_file_to_fastq_file_index($read_file);
 # 	my $this_complex_error_counts_file_name = $settings->file_name('complex_error_counts_file_name', {'#'=> $read_file});
 # 	save_complex_error_file($this_complex_error_counts_file_name, $complex_error_hash->[$fastq_file_index]);
 # }
 
 }
*/