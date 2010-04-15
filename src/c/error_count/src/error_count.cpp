#include <iostream>
#include <cmath>
#include <map>
#include <sstream>
#include <vector>
#include <sam.h>
#include <faidx.h>
#include "error_count.h"


/*! User-defined struct that holds information to be used by the pileup function.
 */
struct user_data {
	
	typedef std::map<std::string,int> base_count_t;
	typedef std::map<uint8_t,base_count_t> qual_map_t;
		
	//! Constructor.
	user_data(const std::string& bam, const std::string& fasta)
	: in(0), ref(0), ref_seq(0), ref_seq_len(0) {
		using namespace std;
		in = samopen(bam.c_str(), "rb", 0);
		ref = fai_load(fasta.c_str());

		// for each sequence id...
		for(int i=0; i<in->header->n_targets; ++i) {
			cerr << "  REFERENCE: " << in->header->target_name[i] << endl;
			cerr << "  LENGTH: " << in->header->target_len[i] << endl;
		}
		
		ostringstream region; region << in->header->target_name[0];
		ref_seq = fai_fetch(ref, region.str().c_str(), &ref_seq_len);
	}

	//! Destructor.
	~user_data() {	
		samclose(in);
		fai_destroy(ref);
		free(ref_seq);
	}

	samfile_t* in; //!< BAM file handle.
	faidx_t* ref; //!< FAI file handle.
	char *ref_seq; //!< Reference sequence.
	int ref_seq_len; //!< Length of the reference sequence.
	
	/*! Coverage count table.
	 
	 This is a table of non-deletion reads per position to non-redundancy counts.
	 For example, given unique_only_coverage[i] = x, for all aligned positions p:
	 i is the number of reads that do not indicate a deletion at p
	 x is the number of positions that have no redundancies
	 */
	std::vector<int> unique_only_coverage;
	
	//my $error_hash = [];		#list by fastq file index
	// if this gets slow, see boost::multi_index_container
	std::map<int32_t,qual_map_t> error_hash; //!< fastq_file_index -> quality map.
};

//1 for A, 2 for C, 4 for G,
//8 for T and 15 for N.
#define is_A(x) (x == 0x01)
#define is_C(x) (x == 0x02)
#define is_G(x) (x == 0x04)
#define is_T(x) (x == 0x08)
#define is_N(x) (x == 0x0f)

/*! Reverse a base.
 */
inline uint8_t reverse_base(uint8_t base) {
	if(base > 0x0f) {
		// ascii
		switch(base) {
			case 'A': return 'T';
			case 'C': return 'G';
			case 'G': return 'C';
			case 'T': return 'A';
			default: assert(false);
		}
	} else {
		// sam-style 4-bit field
		switch(base) {
			case 0x1: return 0x8;
			case 0x2: return 0x4;
			case 0x4: return 0x2;
			case 0x8: return 0x1;
			default: assert(false);
		}		
	}
}

/*! Convert a base to an ASCII character.
 */
inline char base2char(uint8_t base) {
	if(base > 0x0f) {
		// already in ascii format
		return static_cast<char>(base);
	} else {
		// sam-style 4-bit field
		switch(base) {
			case 0x01: return 'A';
			case 0x02: return 'C';
			case 0x04: return 'G';
			case 0x08: return 'T';
			case 0x0f: return '.';
			default: assert(false);
		}		
	}
}


/*! Pileup callback for error count.
 
1) Testing for a relationship between the number of unique reads -> errors.
 
 \param tid target id, corresponding to the index of this target's name
 \param pos position of this alignment in the reference sequence
 \param n number of alignments occurring at this position
 \param pile array of n alignments occuring at this position
 \param data pointer to user_data struct.
 */
int error_count_callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile, void *data) {
	using namespace std;
	user_data* ud = reinterpret_cast<user_data*>(data);
	
	if((pos % 10000) == 0) {
		cerr << "    POSITION:" << pos << endl;
	}

	// reference bases {base, reversed base, null} in ascii:
	char ref_base[] = {ud->ref_seq[pos], reverse_base(ud->ref_seq[pos]), 0};

	size_t unique_coverage=0; // number of non-deletion, non-redundant alignments at this position.
	bool has_redundant_reads=false; // flag indicating whether this position has any redundant reads.
	
	// for each alignment within this pileup:	
	for(int i=0; i<n; ++i) {
		const bam_pileup1_t* p=&pile[i];
		const bam1_t* a=pile[i].b;

		// is this a redundant read?  if so, don't process it - we're all done.
		// also, mark this position as having a redundant read so that we don't update the
		// coverage count when we're done looking at all the alignments.
		bool redundant = (bam_aux2i(bam_aux_get(a,"X1")) > 1);
		if(redundant) {
			has_redundant_reads = true;
			continue;
		}
		
		// track the number of non-deletion, non-redundant alignments:
		if(!p->is_del) {
			++unique_coverage;
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
				base = reverse_base(base);
			}
			
			// my $key = $ref_base[$reversed] . $base; 
			string key; key += static_cast<char>(ref_base[reversed]); key += base2char(base);
			// $error_hash->[$fastq_file_index]->{$quality}->{$key}++;
			++ud->error_hash[fastq_file_index][quality][key];

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
			string key; key += static_cast<char>(ref_base[reversed]); key += '.';
			//	$error_hash->[$fastq_file_index]->{$quality}->{$key}++;	
			++ud->error_hash[fastq_file_index][quality][key];
		} 
		
		if(p->indel == 1) {			
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
			string key; key += '.'; key += base2char(base);
			//	$error_hash->[$fastq_file_index]->{$quality}->{$key}++;
			++ud->error_hash[fastq_file_index][quality][key];
		}
	}
		
	// NOTE: this can't move inside the for-loop; we're tracking information about
	// the position, not each alignment.
	// record coverage at this position, but only if there were no redundant reads:
	if(!has_redundant_reads) {
		// resize our coverage table as needed:
		if(unique_coverage >= ud->unique_only_coverage.size()) {
			ud->unique_only_coverage.resize(unique_coverage+1,0); // >= and +1 because of 0-indexing.
		}
		++ud->unique_only_coverage[unique_coverage];
	}

	return 0;  
}


/*! Print coverage distribution.
 */
void print_coverage_distribution(std::ostream& out, user_data& ud) {
	using namespace std;
	//	# print out coverage, but ZERO OUT observations of zero
	//	# because these may be bona fide deletions
	out << "coverage\tn" << endl;
	//	for (my $i=1; $i<scalar @$unique_only_coverage_list_ref; $i++)
	for(std::size_t i=1; i<ud.unique_only_coverage.size(); ++i) {
		// my $cov = $unique_only_coverage_list_ref->[$i];
		int cov = ud.unique_only_coverage[i];
		// $cov = 0 if (!defined $cov); #some may be undefined, they mean zero
		// print COV "$i\t$cov\n";
		out << i << "\t" << cov << endl;
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
	user_data::qual_map_t& qual_map=ud.error_hash[fastq_file_index];
	
	for(user_data::qual_map_t::reverse_iterator iter=qual_map.rbegin(); iter!=qual_map.rend(); ++iter) {
		out << static_cast<unsigned int>(iter->first);
		for(int i=0; i<5; ++i) {
			for(int j=0; j<5; ++j) {
				string k; k += bases[i]; k += bases[j];
				out << "\t" << iter->second[k];
			}
		}
		out << endl;
	}
}


/*! Calculate the errors in the given BAM file.
 */
void breseq::error_count(const std::string& bam, const std::string& fasta) {
	using namespace std;
	user_data ud(bam, fasta);
	
	sampileup(ud.in, -1, error_count_callback, &ud);

	print_coverage_distribution(cout, ud);
	
	// print_error_file(cout, ud, 0);
}
