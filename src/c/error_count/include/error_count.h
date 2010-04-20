#ifndef _BRESEQ_ERROR_COUNT_H_
#define _BRESEQ_ERROR_COUNT_H_

#include <string>
#include <vector>

namespace breseq {
	
	//! Calculate the errors in the given BAM file.
	void error_count(const std::string& bam, 
									 const std::vector<std::string>& fastas,
									 const std::string& output_dir,
									 const std::string& coverage_suffix,
									 const std::string& error_suffix);
	
} // breseq

#endif
