/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_MUTATION_PREDICTOR_H_
#define _BRESEQ_MUTATION_PREDICTOR_H_

#include "common.h"

#include "alignment.h"
#include "annotated_sequence.h"
#include "genome_diff.h"
#include "settings.h"
#include "fastq.h"
#include "candidate_junctions.h"

using namespace std;

namespace breseq {

	class MutationPredictor
	{
	public:

		static cReferenceSequences ref_seq_info;

		MutationPredictor(cReferenceSequences& ref_seq_info);
		void predict(Settings& settings, genome_diff& gd, uint32_t max_read_length, double avg_read_length = 0.0);

		static bool sort_by_hybrid(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b);
		static bool sort_by_reject_score(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b);
		static bool sort_by_pos(const counted_ptr<diff_entry>& a, const counted_ptr<diff_entry>& b);

	private:

		cSequenceFeature* within_repeat(string seq_id, uint32_t position);

	}; // class MutationPredictor

	inline int32_t n(string input) { return from_string<int32_t>(input); }
	inline bool b(string input) { return from_string(input); }
	inline string s(int32_t input) { return to_string(input); }

} // namespace breseq

#endif
