#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <sam.h>
#include "error_count.h"

typedef std::map<std::string,int> base_count_t;
typedef std::map<uint8_t,base_count_t> qual_map_t;

/*! User-defined struct that holds information to be used by the pileup function.
 */
struct user_data {
	//! Constructor.
	user_data() : in(0) { }
	samfile_t *in; //!< BAM file handle.
	std::vector<int> unique_only_coverage; //!< ?
	// if this gets slow, see boost::multi_index_container
	std::map<int32_t,qual_map_t> error_hash; //!< ?
};

//1 for A, 2 for C, 4 for G,
//8 for T and 15 for N.
#define is_A(x) (x == 0x1)
#define is_C(x) (x == 0x2)
#define is_G(x) (x == 0x4)
#define is_T(x) (x == 0x8)
#define is_N(x) (x == 0xf)

// We're going to play a trick common in FFTs, and just use a lookup table to 
// reverse bases.
uint8_t _reverse_table[9] = {
0x0, // 0:
0x8, // 1: A -> T
0x4, // 2: C -> G
0x0, // 3:
0x2, // 4: G -> C
0x0, // 5:
0x0, // 6:
0x0, // 7:
0x1, // 8: T -> A
};
#define revcom(x) (_reverse_table[x])


/*! Pileup callback for error count.
 
 \param tid chromosome id (?)
 \param pos starting position of this alignment
 \param n number of alignments at this position
 \param pl array of n alignments occuring at this position
 */
int error_count_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
	using namespace std;
	user_data* ud = reinterpret_cast<user_data*>(data);
	
	//my ($seq_id,$pos,$pileup) = @_;
	//print STDERR "    POSITION:$pos\n" if ($pos % 10000 == 0);
	
	if((pos % 10000) == 0) {
		cerr << "    POSITION:" << pos << endl;
	}
	
	//my @ref_base; #index is 'reversed'
	//$ref_base[0] = substr $ref_seq_string, $pos-1, 1;
	//$ref_base[1] = substr $com_seq_string, $pos-1, 1;
	//
	//my $unique_only_position = 1;
	int unique_only_position=1;
	//my $unique_coverage = 0;
	int unique_coverage=0;
	
	// grow our vector of unique coverages as needed, init to 0:
	if(static_cast<std::size_t>(n) > ud->unique_only_coverage.size()) {
		ud->unique_only_coverage.resize(n,0);
	}
	
	//ALIGNMENT: for my $p (@$pileup) 
	for(int i=0; i<n; ++i) {
		const bam_pileup1_t* p=&pl[i];
		
		//my $a = $p->alignment;
		const bam1_t* a=pl[i].b;
		
		//my $indel = $p->indel;
		//$indel = 0 if ($indel < 0);
		//$indel = -1 if ($p->is_del);
		int indel = p->indel;
		if(indel < 0) {
			indel = 0;
		}
		if(p->is_del) {
			indel = -1;
		}
		
		//my $redundancy = $a->aux_get('X1');
		int32_t redundancy=bam_aux2i(bam_aux_get(a,"X1"));
		
		//if (!$p->is_del >= 0)
		//{
		// \todo this should always be true; p->is_del is a bitfield:1...
		
		//	if ($redundancy == 1)
		//	{
		//		$unique_coverage++;
		//	}
		//	else
		//	{
		//		$unique_only_position = 0;
		//	}
		//}
		if(redundancy == 1) {
			++unique_coverage;
		} else {
			unique_only_position = 0;
		}
		
		// record unique only coverage
		// $unique_only_coverage->[$unique_coverage]++ if ($unique_only_position);
		if(unique_only_position) {
			++ud->unique_only_coverage[unique_coverage]; // see above, where we resize this vector.
		}
		
		//next if ($redundancy != 1);
		if(redundancy != 1) {
			continue;	// don't process non-unique reads
		}
		
		//my $qseq = $a->qseq;
		uint8_t* qseq = bam1_seq(a);
		//my $qpos = $p->qpos;
		int32_t qpos = p->qpos;
		//my $qscore = $a->qscore;
		uint8_t* qscore = bam1_qual(a);
		//my $reversed = $a->reversed;
		uint32_t reversed = bam1_strand(a);
		//my $query_start = $a->query->start;
		int32_t query_start = a->core.pos;
		//my $query_end = $a->query->end;
		int32_t query_end = query_start + a->core.l_qseq;
		//my $fastq_file_index = $a->aux_get('X2');
		int32_t fastq_file_index=bam_aux2i(bam_aux_get(a,"X2"));
		
		//In all that follows, be sure to keep track of strandedness of mutations!
		
		if(!p->is_del) {
			//if (!$p->is_del) {
			//# (1) base substitutions
			//#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
			//#     IMPROVE by keeping track of flanking base context and scores?
			//#     this would, for example, penalize low scoring sequences more
			
			// my $base  = substr($qseq, $qpos,1);
			uint8_t base = bam1_seqi(qseq,qpos);
			
			// next if ($base =~ /[nN]/);
			if(is_N(base)) {
				continue;
			}
			
			// my $quality = $qscore->[$qpos];
			uint8_t quality = qscore[qpos];
			
			// $base = FastqLite::revcom($base) if ($reversed);
			if(reversed) {
				base = revcom(base);
			}
			
			// my $key = $ref_base[$reversed] . $base; 
			// $error_hash->[$fastq_file_index]->{$quality}->{$key}++;
			
			// also add an observation of a non-gap non-gap			
			// if ($qpos+1 < $query_end) {
			if((qpos+1) < query_end) {
				// my $next_quality = $qscore->[$qpos+1];
				uint8_t next_quality = qscore[qpos+1];
				
				// my $avg_quality = POSIX::floor( ($quality + $next_quality) / 2);
				uint8_t avg_quality = floor((static_cast<float>(quality)+static_cast<float>(next_quality))/2.0);
				
				// $error_hash->[$fastq_file_index]->{$avg_quality}->{'..'}++;
				++ud->error_hash[fastq_file_index][avg_quality][".."];
			}
		} else if(p->indel == 0) {
			//elsif ($p->indel == 0) {
			//# (2) deletion in read relative to reference
			//#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
			//#     how does one give a quality score? Use quality score of the next base in
			//#     the read, i.e. where the deleted base would have been in the read
			//#     PROBLEM what do we do about multiple base deletions?
			//#     -- only count if there is no indel after this posision
			// train error model only on single base insertions or deletions.
			
			//	my $quality = $qscore->[$qpos+(1-$reversed)];
			uint8_t quality = qscore[qpos+(1-reversed)];
			
			//	my $key = $ref_base[$reversed] . '.';
			
			//	$error_hash->[$fastq_file_index]->{$quality}->{$key}++;	
			
		} else if(p->indel == 1) {			
			//if ($p->indel == +1) {
			//# (3) insertion in read relative to reference
			//#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
			//#     how does one give a quality score? - 
			//#     for reference observations: average the quality scores of the surrounding bases in the read (round down)
			//#     for mutation observations: the quality score of the inserted base
			//#     -- at the next position in the read
			//#     -- only count if an indel = +1, meaning a single-base insertion
			// train error model only on single base insertions or deletions.
			
			//	my $base  = substr($qseq,$qpos+1,1);
			uint8_t base = bam1_seqi(qseq,qpos+1);
			
			//	next if ($base =~ /[nN]/);
			if(is_N(base)) {
				continue;
			}
			
			//	my $quality = $a->qscore->[$qpos+1];
			uint8_t quality = qscore[qpos+1];
			
			//	my $key = '.' . $base; 
			//	$error_hash->[$fastq_file_index]->{$quality}->{$key}++;
		}
	}
	
	return 0;  
}


/*! Print coverage distribution.
 */
void print_coverage_distribution(std::ostream& out, user_data& ud) {
	//	# print out coverage, but ZERO OUT observations of zero
	//	# because these may be bona fide deletions
	//	for (my $i=1; $i<scalar @$unique_only_coverage_list_ref; $i++)
	for(std::size_t i=0; i<ud.unique_only_coverage.size(); ++i) {
		// my $cov = $unique_only_coverage_list_ref->[$i];
		int cov = ud.unique_only_coverage[i];
		// $cov = 0 if (!defined $cov); #some may be undefined, they mean zero
		// print COV "$i\t$cov\n";
		out << i << "\t" << cov << std::endl;
	}	
}


/*! Print error file.
 */
void print_error_file(std::ostream& out, user_data& ud, int32_t fastq_file_index) {
	using namespace std;
	//	#print out a table of errors stratified by quality scores
	//	
	//	##
	//	## Note that we are ignoring 'N' bases by NOT PRINTING THEM OUT
	//	## Later when encountering these they should be SKIPPED
	//	##
	//		
	//	## Create header list which is the same for all files		
	//	my @bases = ('A', 'T', 'C', 'G', '.');
	char bases[] = {'A', 'T', 'C', 'G', '.'};
	//	my @header_list = ('quality');
	out << "quality";
	
	//	foreach my $base_1 (@bases)
	//	{
	//		foreach my $base_2 (@bases)
	//		{
	//			push @header_list, "$base_1$base_2";
	//		}
	//	}
	//	
	//	print ERR join("\t", @header_list). "\n";
	for(int i=0; i<5; ++i) {
		for(int j=0; j<5; ++j) {
			out << "\t" << bases[i] << bases[j];
		}
	}
	out << endl;
	
	//	foreach my $quality (sort { -($a<=>$b) } keys %{$error_hash_ref})
	//	{
	//		my @line_list;
	//		push @line_list, $quality;
	//		
	//		foreach my $base_1 (@bases)
	//		{
	//			foreach my $base_2 (@bases)
	//			{
	//				my $val = $error_hash_ref->{$quality}->{"$base_1$base_2"};
	//				$val = 0 if (!defined $val);
	//				push @line_list, $val;
	//			}
	//		}
	//		print ERR join("\t", @line_list). "\n";
	//	}	
	qual_map_t& qual_map=ud.error_hash[fastq_file_index];
	
	for(qual_map_t::reverse_iterator iter=qual_map.rbegin(); iter!=qual_map.rend(); ++iter) {
		out << static_cast<unsigned int>(iter->first);
		for(int i=0; i<5; ++i) {
			for(int j=0; j<5; ++j) {
				string k; k += bases[i]; k += bases[j];
				out << "\t" << iter->second[k];
			}
		}
		out << endl;
	}
	out << endl;
}


/*! Calculate the errors in the given BAM file.
 */
void breseq::error_count(const std::string& bam) {
	using namespace std;
	
	//my ($settings, $summary, $ref_seq_info) = @_;
	//
	//my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');
	//my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	//my $bam = Bio::DB::Sam->new(-fasta => $reference_fasta_file_name, -bam => $reference_bam_file_name);
	//my @seq_ids = $bam->seq_ids;
	//
	//## populated by pileup
	//my $error_hash = [];		#list by fastq file index
	
	user_data ud;
	ud.in = samopen(bam.c_str(), "rb", 0);
	
	//foreach my $seq_id (@seq_ids)
	//{							
	//	my $sequence_length = $bam->length($seq_id);
	//	print STDERR "  REFERENCE: $seq_id\n";
	//	print STDERR "  LENGTH: $sequence_length\n";
	
	// for each sequence id...
	for(int i=0; i<ud.in->header->n_targets; ++i) {
		cerr << "  REFERENCE: " << ud.in->header->target_name[i] << endl;
		cerr << "  LENGTH: " << ud.in->header->target_len[i] << endl;
	}
	
	//## populated by pileup
	//my $unique_only_coverage;
	//my $ref_seq_string = $ref_seq_info->{ref_strings}->{$seq_id};
	//my $com_seq_string = $ref_seq_string;
	//$com_seq_string =~ tr/ATCG/TAGC/;
	
	// $bam->pileup($seq_id,$pileup_function);
	sampileup(ud.in, -1, error_count_callback, &ud);
	samclose(ud.in);
	
	// print_coverage_distribution(cout, ud);
	
	print_error_file(cout, ud, 0);
}
