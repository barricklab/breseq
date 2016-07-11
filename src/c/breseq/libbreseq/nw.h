
/*
 * nw.h for program nw.
 *
 * Modified for breseq from http://www.rolfmuertter.com/code/nw.php
 *
 */

#ifndef _BRESEQ_NW_H_
#define _BRESEQ_NW_H_

#include "common.h"

#define DEBUG 0

using namespace std;

namespace breseq {


extern int  nw( 
                 string, string, 
                 string&, string&,
                 bool, int32_t
              );

extern int32_t  nw_align(
                      int32_t **, char **,
                      string, string, 
                      string&, string&,
                      int32_t, int32_t
                    );

extern void  dpm_init        ( int32_t **, char **, int32_t, int32_t, int32_t );
extern void  print_al        ( string&, string& );
extern void  print_matrix    ( int32_t ** const, string, string );
extern void  print_traceback ( char ** const, string, string );
extern int32_t   nw_max             ( int32_t, int32_t, int32_t, char * );

} // namespace breseq

#endif

