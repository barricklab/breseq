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
#include "breseq/pileup.h"
#include "breseq/error_count.h"


/*! Count errors.
 */
void breseq::error_count(const std::string& bam, 
												 const std::string& fasta,
												 const std::string& output_dir,
												 const std::vector<std::string>& readfiles,
                         bool do_coverage,
                         bool do_errors,
                         uint8_t min_qual_score
                        ) 
{
	error_count_pileup ecp(bam, fasta, do_coverage, do_errors, min_qual_score);
	ecp.do_pileup();
	if (do_coverage) ecp.print_coverage(output_dir);
	if (do_errors) ecp.print_error(output_dir, readfiles);
}


/*! Constructor.
 */
breseq::error_count_pileup::error_count_pileup(const std::string& bam, const std::string& fasta, bool do_coverage, bool do_errors, uint8_t min_qual_score)
: breseq::pileup_base(bam, fasta), m_do_coverage(do_coverage), m_do_errors(do_errors), m_min_qual_score(min_qual_score) {
	// reserve enough space for the sequence info:
	_seq_info.resize(_bam->header->n_targets);
}


/*! Destructor.
 */
breseq::error_count_pileup::~error_count_pileup() {
}


/*! Called for each alignment.
 */
void breseq::error_count_pileup::callback(const breseq::pileup& p) {
	using namespace std;
	assert(p.target() < _seq_info.size());
	sequence_info& info=_seq_info[p.target()];

	size_t unique_coverage=0; // number of non-deletion, non-redundant alignments at this position.
	bool has_redundant_reads=false; // flag indicating whether this position has any redundant reads.

	// for each alignment within this pileup:
	for(pileup::const_iterator i=p.begin(); i!=p.end(); ++i) {
		
    // skip deletions entirely, they are handled by adjacent matching positions
    if(i->is_del()) {
      continue;
    }
		
		// is this a redundant read?  if so, don't process it - we're all done.
		// also, mark this position as having a redundant read so that we don't update the
		// coverage count when we're done looking at all the alignments.
		if(i->is_redundant()) {
			has_redundant_reads = true;
			continue;
		}
		
		// track the number of non-deletion, non-redundant alignments:
		++unique_coverage;
		
		// if we are only tracking coverage, go to the next alignment
    if (!m_do_errors) {
      continue;
    }		
				
		uint32_t reversed = i->strand(); // are we on the reverse strand?
		uint8_t* qseq = i->query_sequence(); // query sequence (read)
		int32_t qpos = i->query_position(); // position of the alignment in the query

		int32_t qstart = i->query_start() - 1; // @dk: 0-indexed, so subtract 1??  (should add)
		int32_t qend = i->query_end() - 1;

		int32_t qlen = i->query_length(); // length of this query
		uint8_t* qscore = i->quality_scores(); // quality score array
		int32_t fastq_file_index = i->fastq_file_index(); // sequencer-generated read file that this alignment belongs to
		
		uint32_t pos = p.position(); // position of this alignment on the reference sequence
		char* refseq = p.reference_sequence(); // reference sequence for this target
		char ref_base[] = {refseq[pos], reverse_base(refseq[pos]), 0}; // reference base & its complement
		
		// Things to remember in the following:
    // -->1 Reverse the base when the read is on the other strand
    // -->2 Correct for which base is AFTER in error counts when on the other strand
    // -->3 Observations that involve an N base in any way should not be counted

			//# (1) base substitutions
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
        ++_error_hash[fastq_file_index][quality][key];
      }
    }

		//# the next base also matches 
    //# (1) base substitution or match
    //#     e.g. '..' key indicating an observation of a "non-gap, non-gap"
    //#     quality score is of the second non-gap in the pair
    if(i->indel() == 0) {
      //## don't count past last match position
      if (qpos < qend) {	
        int32_t mqpos = qpos + 1 - reversed;
        uint8_t base = bam1_seqi(qseq,mqpos);
        
        int32_t mrpos = pos + 1 - reversed;
        uint8_t ref_base = refseq[mrpos];
        
        if (!is_N(base) && !is_char_N(ref_base)) {     
          string key = ".."; 
          uint8_t quality = qscore[mqpos];
          ++_error_hash[fastq_file_index][quality][key];
        }
      }	
    }
		
		//# there is a deletion of EXACTLY one base in the read relative to the reference before the next read base
    //# (2) deletion in read relative to reference
    //#     e.g. 'A.' key for observing nothing in a read at a position where the reference has an A
    //#     quality score is of the next non-gap base in the read
    else if (i->indel() == -1) {
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
        ++_error_hash[fastq_file_index][quality][key];									
      }
    }

		
    //# there is an insertion of EXACTLY one base in the read relative to the reference before the next reference base
    //# (3) insertion in read relative to reference
    //#     e.g. '.A' key for observing an A in a read at a position where the reference has no base
    //#     quality score is that of the observed inserted base
    else if (i->indel() == +1) {
      int32_t mqpos = qpos + 1;
			
      if ((mqpos <= qend) && (mqpos >= qstart)) {
        uint8_t base = bam1_seqi(qseq,mqpos);    
				
        if (!is_N(base)) {
          if(reversed) {
            base = reverse_base(base);
          }        
          string key; key += '.'; key += base2char(base);
          uint8_t quality = qscore[mqpos];
          ++_error_hash[fastq_file_index][quality][key];									
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
	
	assert(readfiles.size() == _error_hash.size());
	
	for(fastq_map_t::iterator iter=_error_hash.begin(); iter!=_error_hash.end(); ++iter) {
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
      if (iter->first < m_min_qual_score) continue; // skip low quality scores
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


/*! Load error rates.
 */
breseq::error_count_results::error_count_results(const std::string& input_dir, const std::vector<std::string>& readfiles) {
	using namespace std;
	char bases[] = {'A', 'T', 'C', 'G', '.'}; // order is important!!! must match the header from the error rates file...
	
	_error_rates.resize(readfiles.size());
	
	// load the error rates:
	for(std::size_t i=0; i<readfiles.size(); ++i) {
		string filename(input_dir + readfiles[i] + ".error_rates.tab");
		ifstream in(filename.c_str());
		in.ignore(1024, '\n'); // get rid of header line.
		
		error_map_t& emt = _error_rates[i];
		
		while(!in.eof()) {
			int quality;
			in >> quality;
			
			error_map_t::iterator em = emt.insert(make_pair(static_cast<uint8_t>(quality),base_error_t())).first;
			
			for(int i=0; i<5; ++i) {
				for(int j=0; j<5; ++j) {
					string k; k += bases[i]; k += bases[j];
					double r;
					in >> r;
					double cor = (r==0.0 ? 1e-256 : log10(r));
					double err = (r==1.0 ? 1e-256 : log10(1.0-r));
					em->second[k] = make_pair(err,cor);
				}
			}
		}
		in.close();
	}
}

/*! Return the error rate for the given base pair, quality, and FASTQ file index.
 */
double breseq::error_count_results::log10_error_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
	//  fail if this doesn't exist?
	assert((std::size_t)fastq_file_index<_error_rates.size());
	assert(_error_rates[fastq_file_index].find(quality)!=_error_rates[fastq_file_index].end());
	assert(_error_rates[fastq_file_index][quality].find(base_key)!=_error_rates[fastq_file_index][quality].end());
	
	return _error_rates[fastq_file_index][quality][base_key].first;
}


/*! Return the correct rate for the given base pair, quality, and FASTQ file index.
 */
double breseq::error_count_results::log10_correct_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
		//  fail if this doesn't exist?
	assert((std::size_t)fastq_file_index<_error_rates.size());
	assert(_error_rates[fastq_file_index].find(quality)!=_error_rates[fastq_file_index].end());
	assert(_error_rates[fastq_file_index][quality].find(base_key)!=_error_rates[fastq_file_index][quality].end());
	
	return _error_rates[fastq_file_index][quality][base_key].second;
}

/*! Return the pair of (correct,error) rates.
 */
const std::pair<double,double>& breseq::error_count_results::log10_rates(int32_t fastq_file_index, uint8_t quality, const std::string& base_key) {
	return _error_rates[fastq_file_index][quality][base_key];
}
