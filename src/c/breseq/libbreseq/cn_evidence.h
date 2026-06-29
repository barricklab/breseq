/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2026 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _BRESEQ_CN_EVIDENCE_H_
#define _BRESEQ_CN_EVIDENCE_H_

#include "common.h"

#include "reference_sequence.h"
#include "settings.h"

namespace breseq {

  // Runs the external tool CNery (https://github.com/barricklab/CNery) on the
  // current breseq output, then ingests its per-reference-sequence copy
  // number calls into the CN evidence genome diff files that the rest of the
  // pipeline (mutation prediction passthrough, HTML display) already expects.
  class CNEvidence
  {
  public:

    static void predict(
                        Settings& settings,
                        Summary& summary,
                        cReferenceSequences& ref_seq_info
                        );

  private:

    static void run_cnery(Settings& settings, Summary& summary, const string& cnery_output_prefix);

    static void ingest_csv_for_seq_id(
                                      const string& seq_id,
                                      const string& cnv_file_name,
                                      const string& break_pts_file_name,
                                      const string& gd_file_name
                                      );

  }; // class CNEvidence

} // namespace breseq

#endif
