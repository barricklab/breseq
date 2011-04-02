/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2010 Michigan State University

  breseq is free software; you can redistribute it and/or modify it under the  
  terms the GNU General Public License as published by the Free Software 
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#ifndef _CANDIDATE_JUNCTIONS_H_
#define _CANDIDATE_JUNCTIONS_H_

#include <assert.h>
#include <string>

#include <bam.h>

#include "common.h"

using namespace std;

namespace breseq {

  /*! Predicts candidate junctions
   */
  void candidate_junctions(string fasta, string sam, string output);
	
} // namespace breseq

#endif
