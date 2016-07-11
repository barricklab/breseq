/*-------------------------------------------
 * 
 *  nw.c++ for program nw
 *
 *  Modified for breseq from http://www.rolfmuertter.com/code/nw.php
 * 
 -------------------------------------------*/

#include "libbreseq/nw.h"

using namespace std;

namespace breseq {

int nw(                                                          
        string       seq_1,          /*  Needleman-Wunsch   */
        string       seq_2,          /*  algorithm for      */
        string&      seq_1_al,       /*  global alignment   */
        string&      seq_2_al,       /*  of nt sequence.    */
        bool         prm,
        int32_t      min_score
      )
{
  int  d = -1 ;                 /* gap score */

  int  L1 = seq_1.length();
  int  L2 = seq_2.length();

  // Dynamic programming matrix
  int ** F = new int32_t * [ L2+1 ];
  for( int32_t i = 0; i <= L2; i++ )  F[ i ] = new int [ L1+1 ];

  // Traceback matrix
  char ** traceback = new char * [ L2+1 ];
  for( int32_t i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1+1 ];

  // Initialize traceback and F matrix (fill in first row and column)
  dpm_init( F, traceback, L1, L2, d );

  // Create alignment
  int32_t score = nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d, min_score );

  #if DEBUG
      int  L_al = seq_1_al.length();
      cout << "Length after alignment: " << L_al << endl;
  #endif

  if( prm )
  {
    cout << "\nDynamic programming matrix: " << "\n\n";
    print_matrix( F, seq_1, seq_2 );

    cout << "\nTraceback matrix: " << "\n\n";
    print_traceback( traceback, seq_1, seq_2 );

    cout << endl;
  }

  for( int i = 0; i <= L2; i++ )  delete F[ i ];
  delete[] F;
  for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
  delete[] traceback;

  return  score ;
}


void  dpm_init( int32_t ** F, char ** traceback, int32_t L1, int32_t L2, int32_t d )
{
  F[ 0 ][ 0 ] =  0 ;
  traceback[ 0 ][ 0 ] = 'n' ;

  int i=0, j=0;

  for( j = 1; j <= L1; j++ ) {
    F[ 0 ][ j ] =  j * d ;
    traceback[ 0 ][ j ] =  '-' ;
  }
  
  for( i = 1; i <= L2; i++ ) {
    F[ i ][ 0 ] =  i * d ;
    traceback[ i ][ 0 ] =  '|' ;
  }
}


int nw_align(                  // Needleman-Wunsch algorithm
              int32_t **     F,
              char **    traceback,
              string     seq_1,
              string     seq_2,
              string&    seq_1_al,
              string&    seq_2_al,
              int32_t    d,
              int32_t    min_score
            )
{
  seq_1_al = "";
  seq_2_al = "";
  
  
  int32_t        k = 0, x = 0, y = 0;
  int32_t        fU, fD, fL ;
  char       ptr, nuc ;
  int32_t        i = 0, j = 0;

  const int32_t  a =  0;   // Match
  const int32_t  b = -1;   // Mismatch


  const int32_t  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                             { b, a, b, b },
                             { b, b, a, b },
                             { b, b, b, a } } ;

  int32_t  L1 = seq_1.length();
  int32_t  L2 = seq_2.length();

  for( i = 1; i <= L2; i++ )
  {
    int32_t max_row_score = numeric_limits<int32_t>::min();
    for( j = 1; j <= L1; j++ )
    {
      nuc = seq_1[ j-1 ] ;

      switch( nuc )
      {
        case 'A':  x = 0 ;  break ;
        case 'C':  x = 1 ;  break ;
        case 'G':  x = 2 ;  break ;
        case 'T':  x = 3 ;
      }

      nuc = seq_2[ i-1 ] ;

      switch( nuc )
      {
        case 'A':  y = 0 ;  break ;
        case 'C':  y = 1 ;  break ;
        case 'G':  y = 2 ;  break ;
        case 'T':  y = 3 ;
      }

      fU = F[ i-1 ][ j ] + d ;
      fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
      fL = F[ i ][ j-1 ] + d ;

      F[ i ][ j ] = nw_max( fU, fD, fL, &ptr ) ;
      
      max_row_score = max(max_row_score, F[ i ][ j ]);

      traceback[ i ][ j ] =  ptr ;
    }
    
    // Leave early if we are above the max, this is not as efficient as would be possible
    if (max_row_score < min_score) return max_row_score;
  }
  i-- ; j-- ;
  int32_t score = F[ i ][ j ];

  while( i > 0 || j > 0 )
  {
    switch( traceback[ i ][ j ] )
    {
      case '|' :      seq_1_al += '-' ;
                      seq_2_al += seq_2[ i-1 ] ; 
                      i-- ;
                      break ;

      case '\\':      seq_1_al += seq_1[ j-1 ] ; 
                      seq_2_al += seq_2[ i-1 ] ; 
                      i-- ;  j-- ;
                      break ;

      case '-' :      seq_1_al += seq_1[ j-1 ] ; 
                      seq_2_al += '-' ; 
                      j-- ;
    }
    k++ ;
  }
  
  reverse( seq_1_al.begin(), seq_1_al.end() );
  reverse( seq_2_al.begin(), seq_2_al.end() );

  return  score ;
}


int32_t  nw_max( int32_t f1, int32_t f2, int f3, char * ptr )
{
  int32_t  max = 0 ;

  if( f1 >= f2 && f1 >= f3 ) {
    max = f1 ;
    *ptr = '|' ;
  } else if( f2 > f3 ) {
    max = f2 ;
    *ptr = '\\' ;
  } else {
    max = f3 ;
    *ptr = '-' ;
  }
  
  return  max ;
}


void  print_matrix( int32_t ** F, string seq_1, string seq_2 )
{
  int32_t  L1 = seq_1.length();
  int32_t  L2 = seq_2.length();

  cout << "        ";
  for( int32_t j = 0; j < L1; j++ ) {
    cout << seq_1[ j ] << "   ";
  }
  cout << "\n  ";

  for( int32_t i = 0; i <= L2; i++ ) {
    
    if( i > 0 ) {
      cout << seq_2[ i-1 ] << " ";
    }
    
    for( int32_t j = 0; j <= L1; j++ ) {
      cout.width( 3 );
      cout << F[ i ][ j ] << " ";
    }
    
    cout << endl;
  }
}


void  print_traceback( char ** traceback, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "    ";
        for( int32_t j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";

        for( int32_t i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int32_t j = 0; j <= L1; j++ )
                {
                        cout << traceback[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  print_al( string& seq_1_al, string& seq_2_al )
{
  cout << seq_1_al << endl;
  cout << seq_2_al << endl;
}

} // namespace breseq
