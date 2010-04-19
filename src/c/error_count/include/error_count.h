#ifndef _BRESEQ_ERROR_COUNT_H_
#define _BRESEQ_ERROR_COUNT_H_

#include <string>

namespace breseq {
	
	//! Calculate the errors in the given BAM file.
	void error_count(const std::string& bam, const std::string& fasta, const std::string& output);
	
} // breseq

#endif
