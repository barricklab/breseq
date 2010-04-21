#include <iostream>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <sam.h>

#include "common.h"
#include "error_count.h"
#include "pileup.h"


namespace breseq {
	
	/*! Error-counting class.
	 */
	class error_count_pileup : public breseq::pileup_base {
	public:
		typedef std::map<std::string,int> base_count_t;
		typedef std::map<uint8_t,base_count_t> qual_map_t;
		
		/*! Information that is tracked per-sequence.
		 */
		struct sequence_info {
			/*! Coverage count table.
			 
			 This is a table of non-deletion reads per position to non-redundancy counts.
			 For example, given unique_only_coverage[i] = x, for all aligned positions p:
			 i is the number of reads that do not indicate a deletion at p
			 x is the number of positions that have no redundancies
			 */
			std::vector<int> unique_only_coverage;
			
			/*! fastq_file_index -> quality map.
			 */
			std::map<int32_t,qual_map_t> error_hash;
		};
		
		
		/*! Constructor.
		 */
		error_count_pileup(const std::string& bam, const std::vector<std::string>& fastas)
		: breseq::pileup_base(bam, fastas) {
			// reserve enough space for the sequence info:
			_seq_info.resize(_bam->header->n_targets);
		}
		
		
		/*! Destructor.
		 */
		virtual ~error_count_pileup() {
		}
		
		
		/*! Called for each alignment.
		 */
		virtual int callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile) {
			using namespace std;
			assert(static_cast<std::size_t>(tid) < _seq_info.size());
			sequence_info& info=_seq_info[tid];
			
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
				
				char* refseq = get_refseq(tid, fastq_file_index);																	
				char ref_base[] = {refseq[pos], reverse_base(refseq[pos]), 0};
				
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
					++info.error_hash[fastq_file_index][quality][key];
					// also add an observation of a non-gap non-gap			
					// if ($qpos+1 < $query_end) {
					if((qpos+1) < query_end) {
						// my $next_quality = $qscore->[$qpos+1];
						uint8_t next_quality = qscore[qpos+1];
						
						// my $avg_quality = POSIX::floor( ($quality + $next_quality) / 2);
						uint8_t avg_quality = floor((static_cast<float>(quality)+static_cast<float>(next_quality))/2.0);
						
						// $error_hash->[$fastq_file_index]->{$avg_quality}->{'..'}++;
						++info.error_hash[fastq_file_index][avg_quality][".."];
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
					++info.error_hash[fastq_file_index][quality][key];
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
					++info.error_hash[fastq_file_index][quality][key];
				}
			}
			
			// NOTE: this can't move inside the for-loop; we're tracking information about
			// the position, not each alignment.
			// record coverage at this position, but only if there were no redundant reads:
			if(!has_redundant_reads) {
				// resize our coverage table as needed:
				if(unique_coverage >= info.unique_only_coverage.size()) {
					info.unique_only_coverage.resize(unique_coverage+1,0); // >= and +1 because of 0-indexing.
				}
				++info.unique_only_coverage[unique_coverage];
			}
			
			return 0;  
		}
		
		
		/*! Print coverage distribution.
		 */
		void print_coverage(const std::string& output_dir, const std::string& suffix) {
			using namespace std;
			for(std::size_t i=0; i<_seq_info.size(); ++i) {
				string filename(output_dir + _bam->header->target_name[i] + "." + suffix);
				ofstream out(filename.c_str());					
				out << "coverage\tn" << endl;
				for(std::size_t j=1; j<_seq_info[i].unique_only_coverage.size(); ++j) {
					out << j << "\t" << _seq_info[i].unique_only_coverage[j] << endl;
				}	
				out.close();
			}
		}
		
		
		/*! Print error file.
		 */
		void print_error(const std::string& output_dir, const std::string& suffix) {
			//										 
			//										 //std::ostream& out, int32_t fastq_file_index) {
			//			using namespace std;
			//			char bases[] = {'A', 'T', 'C', 'G', '.'};
			//			out << "quality";
			//			
			//			for(int i=0; i<5; ++i) {
			//				for(int j=0; j<5; ++j) {
			//					out << "\t" << bases[i] << bases[j];
			//				}
			//			}
			//			out << endl;
			//			
			//			qual_map_t& qual_map=error_hash[fastq_file_index];			
			//			for(qual_map_t::reverse_iterator iter=qual_map.rbegin(); iter!=qual_map.rend(); ++iter) {
			//				out << static_cast<unsigned int>(iter->first);
			//				for(int i=0; i<5; ++i) {
			//					for(int j=0; j<5; ++j) {
			//						string k; k += bases[i]; k += bases[j];
			//						out << "\t" << iter->second[k];
			//					}
			//				}
			//				out << endl;
			//			}
		}
		
		//		for(std::size_t i=0; i<fasta.size(); ++i) {
		//			
		//			if(output.empty()) {
		//				ecp.print_coverage_distribution(cout);
		//				ecp.print_error_file(cout, 0);
		//			} else {
		//				string filename(output+"/unique_only_coverage_distribution.tab");
		//				ofstream coverage(filename.c_str());					
		//				ecp.print_coverage_distribution(coverage);
		//				coverage.close();
		//				
		//				filename = output+"/error_counts.tab";
		//				ofstream errors(filename.c_str());
		//				ecp.print_error_file(errors, 0);
		//				errors.close();
		//			}
		//		}
		// for each sequence id (target id) in the bam file, print the coverage distribution
		// for each reference sequence, print the error counts
		
	protected:		
		std::vector<sequence_info> _seq_info;
	};
	
	
	void error_count(const std::string& bam, 
									 const std::vector<std::string>& fastas,
									 const std::string& output_dir,
									 const std::string& coverage_suffix,
									 const std::string& error_suffix) {
		error_count_pileup ecp(bam, fastas);
		ecp.pileup();
		ecp.print_coverage(output_dir, coverage_suffix);
		ecp.print_error(output_dir, error_suffix);
	}
	
} // breseq
