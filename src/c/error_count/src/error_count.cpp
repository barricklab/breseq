#include <iostream>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <sam.h>
#include <faidx.h>
#include <assert.h>

#include "breseq/common.h"
#include "error_count.h"


/*! Constructor.
 */
breseq::error_count_pileup::error_count_pileup(const std::string& bam, const std::string& fasta, bool do_coverage, bool do_errors)
: breseq::pileup_base(bam, fasta), m_do_coverage(do_coverage), m_do_errors(do_errors) {
	// reserve enough space for the sequence info:
	_seq_info.resize(_bam->header->n_targets);
}


/*! Destructor.
 */
breseq::error_count_pileup::~error_count_pileup() {
}


/*! Called for each alignment.
 */
int breseq::error_count_pileup::callback(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pile) {
	using namespace std;
	assert(static_cast<std::size_t>(tid) < _seq_info.size());
	sequence_info& info=_seq_info[tid];
	
	size_t unique_coverage=0; // number of non-deletion, non-redundant alignments at this position.
	bool has_redundant_reads=false; // flag indicating whether this position has any redundant reads.
	
	// for each alignment within this pileup:	
	for(int i=0; i<n; ++i) {
		const bam_pileup1_t* p=&pile[i];
		const bam1_t* a=pile[i].b;
    
    // skip deletions entirely, they are handled by adjacent matching positions
    if(p->is_del) {
      continue;
    }
		
		// is this a redundant read?  if so, don't process it - we're all done.
		// also, mark this position as having a redundant read so that we don't update the
		// coverage count when we're done looking at all the alignments.
		bool redundant = (bam_aux2i(bam_aux_get(a,"X1")) > 1);
		if(redundant) {
			has_redundant_reads = true;
			continue;
		}
		
		// track the number of non-deletion, non-redundant alignments to this reference position:
    ++unique_coverage;
		
    // if we are only tracking coverage, go to the next alignment
    if (!m_do_errors) {
      continue;
    }
    
		uint32_t reversed = bam1_strand(a); // are we on the reverse strand?
		uint8_t* qseq = bam1_seq(a); // query sequence (read)
		int32_t qpos = p->qpos; // position of the alignment in the query

		int32_t qstart = query_start(a) - 1; //0-indexed
		int32_t qend = query_end(a) - 1; //0-indexed  

		uint8_t* qscore = bam1_qual(a); // quality score array
		int32_t fastq_file_index=bam_aux2i(bam_aux_get(a,"X2")); // sequencer-generated read file that this alignment belongs to
		
		char* refseq = get_refseq(tid); // reference sequence for this target

    // Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

    //# (1) base substitution or match
    //#     e.g. 'AG' key for observing a G in read at a place where the reference has an A
    //#     IMPROVE by keeping track of flanking base context and scores?
    //#     this would, for example, penalize low scoring sequences more
    {
      uint8_t base = bam1_seqi(qseq,qpos);    
      uint8_t ref_base = refseq[pos];
      
      if (!is_N(base) && !is_char_N(ref_base)) {
        if(reversed) {
          base = reverse_base(base);
          ref_base = reverse_base(ref_base);
        }
        string key; key += static_cast<char>(ref_base); key += base2char(base);
        uint8_t quality = qscore[qpos];
        ++error_hash[fastq_file_index][quality][key];
      }
    }
    
    //# the next base also matches 
    //# (1) base substitution or match
    //#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
    //#     quality score is of the second non-gap in the pair
    if (p->indel == 0) {
      //## don't count past last match position
      if (qpos < qend) {	
        int32_t mqpos = qpos + 1 - reversed;
        uint8_t base = bam1_seqi(qseq,mqpos);
        
        int32_t mrpos = pos + 1 - reversed;
        uint8_t ref_base = refseq[mrpos];
        
        if (!is_N(base) && !is_char_N(ref_base)) {     
          string key = ".."; 
          uint8_t quality = qscore[mqpos];
          ++error_hash[fastq_file_index][quality][key];
        }
      }	
    }
	
    //# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
    //# (2) deletion in read relative to reference
    //#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
    //#     quality score is of the next non-gap base in the read
    else if (p->indel == -1) {
      //## count the quality of this or next base depending on reversed, and make sure it is not an N
      int32_t mqpos = qpos + 1 - reversed;
      uint8_t base = bam1_seqi(qseq,mqpos);   
         
      //## the reference base opposite the deletion is really the NEXT base
      int32_t mrpos = pos + 1;
      uint8_t ref_base = refseq[mrpos];      
      
      if (!is_N(base) && !is_char_N(ref_base)) {
        if(reversed) {
          ref_base = reverse_base(ref_base);
        }
        string key; key += static_cast<char>(ref_base); key += '.';
        uint8_t quality = qscore[mqpos];
        ++error_hash[fastq_file_index][quality][key];									
      }
    }
    
    //# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
    //# (3) insertion in read relative to reference
    //#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
    //#     quality score is that of the observed inserted base
    else if (p->indel == +1) {
      int32_t mqpos = qpos + 1;
    
      if ((mqpos <= qend) && (mqpos >= qstart)) {
        uint8_t base = bam1_seqi(qseq,mqpos);    

        if (!is_N(base)) {
          if(reversed) {
            base = reverse_base(base);
          }        
          string key; key += '.'; key += base2char(base);
          uint8_t quality = qscore[mqpos];
          ++error_hash[fastq_file_index][quality][key];									
        }	
      }
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
void breseq::error_count_pileup::print_coverage(const std::string& output_dir) {
	using namespace std;
	for(std::size_t i=0; i<_seq_info.size(); ++i) {
		string filename(output_dir + _bam->header->target_name[i] + ".unique_only_coverage_distribution.tab");
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
void breseq::error_count_pileup::print_error(const std::string& output_dir, const std::vector<std::string>& readfiles) {
	using namespace std;
	char bases[] = {'A', 'T', 'C', 'G', '.'};
	
	assert(readfiles.size() == error_hash.size());
	
	for(fastq_map_t::iterator iter=error_hash.begin(); iter!=error_hash.end(); ++iter) {
		ostringstream filename;
		filename << output_dir << readfiles[iter->first] << ".error_counts.tab";
		ofstream out(filename.str().c_str());
		
		out << "quality";
		for(int i=0; i<5; ++i) {
			for(int j=0; j<5; ++j) {
				out << "\t" << bases[i] << bases[j];
			}
		}
		out << endl;
		
		qual_map_t& qual_map=iter->second;
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
		out.close();
	}
}


/*!
 */
void breseq::error_count(const std::string& bam, 
												 const std::string& fasta,
												 const std::string& output_dir,
												 const std::vector<std::string>& readfiles,
                         const bool do_coverage,
                         const bool do_errors) {
	error_count_pileup ecp(bam, fasta, do_coverage, do_errors);
	ecp.pileup();
	if (do_coverage) ecp.print_coverage(output_dir);
	if (do_errors) ecp.print_error(output_dir, readfiles);
}
