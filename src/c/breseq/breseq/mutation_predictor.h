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

using namespace std;

namespace breseq {

	class MutationPredictor
	{
	public:

		static cReferenceSequences ref_seq_info;

		MutationPredictor(cReferenceSequences& ref_seq_info);
		void predict(Settings& settings, const Summary& summary, genome_diff& gd);

		static bool sort_by_hybrid(diff_entry a, diff_entry b);
		static bool sort_by_reject_score(diff_entry a, diff_entry b);
		static bool sort_by_pos(diff_entry a, diff_entry b);

	private:

		string get_sequence(string seq_id, int32_t start, int32_t end);
		cSequenceFeature* within_repeat(string seq_id, uint32_t position);

	}; // class MutationPredictor

	int32_t n(string input) { return from_string<int32_t>(input); };
	bool b(string input) { return from_string<bool>(input); };
	string s(int32_t input) { return to_string(input); };

} // namespace breseq

#endif
