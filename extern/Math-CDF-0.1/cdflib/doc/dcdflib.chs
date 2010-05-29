










                                    DCDFLIB

               Library of C Routines for Cumulative Distribution
                 Functions, Inverses, and Other Parameters

                                  Version 1.1

                                (November, 1997)






                    Summary Documentation of Each Routine






                            Compiled and Written by:

                                 Barry W. Brown
                                  James Lovato
                                  Kathy Russell









                     Department of Biomathematics, Box 237
                     The University of Texas, M.D. Anderson Cancer Center
                     1515 Holcombe Boulevard
                     Houston, TX      77030


 This work was supported by grant CA-16672 from the National Cancer Institute.


WHICH  and  STATUS are  pointers  to  int ,  all other arguements are
pointers to double.


-------------------- DISTRIBUTION                 

WHICH   PARAMETERS     INPUT RANGE       SEARCH RANGE       REQUIREMENTS


-------------------- Beta 

void cdfbet(int *which,double *p,double *q,double *x,double *y,
            double *a,double *b,int *status,double *bound)

  1     P and Q        [0,1];[0,1]       -----------        SUM to 1.0
  2     X and Y        [0,1];[0,1]       [0,1],[0,1]        SUM to 1.0
  3     A              (0,infinity)      [1E-100,1E100]
  4     B              (0,infinity)      [1E-100,1E100]

-------------------- Binomial

void cdfbin(int *which,double *p,double *q,double *s,double *xn,
            double *pr,double *ompr,int *status,double *bound)

  1     P and Q        [0,1];[0,1]       -----------        SUM to 1.0
  2     S              [0,XN]            [0,XN]
  3     XN             (0,infinity)      [1E-100,1E100]
  4     PR and OMPR    [0,1];[0,1]       [0,1];[0,1]        SUM to 1.0

-------------------- Chi-square

void cdfchi(int *which,double *p,double *q,double *x,double *df,
            int *status,double *bound)

      SUBROUTINE CDFCHI( WHICH, P, Q, X, DF, STATUS, BOUND )

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     X              [0,infinity]      [0,1E100]
  3     DF             (0,infinity)      [1E-100,1E100]

-------------------- Noncentral Chi-square

void cdfchn(int *which,double *p,double *q,double *x,double *df,
            double *pnonc,int *status,double *bound)

  1     P and Q        [0,1-1E-16],none  -----------
  2     X              [0,infinity]      [0,1E100]
  3     DF             (0,infinity)      [1E-100,1E100]
  4     PNONC          [0,infinity)      [0,1E4]

NOTE: We do not yet have a method to calculation the Noncentral Chi-Square
distribution acurately near 1;  therefore, Q is not used by CDFCHN.  There
are no input requirements of Q, and when WHICH is 1, Q is returned as 1-P.

-------------------- F

void cdff(int *which,double *p,double *q,double *f,double *dfn,
          double *dfd,int *status,double *bound)

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     F              [0,infinity)      [0,1E100]
  3     DFN            (0,infinity)      [1E-100,1E100]
  4     DFD            (0,infinity)      [1E-100,1E100]

-------------------- Noncentral F

void cdffnc(int *which,double *p,double *q,double *f,double *dfn,
            double *dfd,double *phonc,int *status,double *bound)

  1     P and Q        [0,1-1E-16],none  -----------
  2     F              [0,infinity)      [0,1E100]
  3     DFN            (0,infinity)      [1E-100,1E100]
  4     DFD            (0,infinity)      [1E-100,1E100]
  5     PNONC          [0,infinity)      [0,1E4]

NOTE: We do not yet have a method to calculation the Noncentral F
distribution acurately near 1;  therefore, Q is not used by CDFF. 
There are no input requirements of Q, and when WHICH is 1, Q is returned
as 1-P.

-------------------- Gamma

void cdfgam(int *which,double *p,double *q,double *x,
            double *shape,double *scale,int *status,double *bound)

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     X              [0,infinity)      [0,1E100]
  3     SHAPE          (0,infinity)      [1E-100,1E100]
  4     SCALE          (0,infinity)      [1E-100,1E100]

-------------------- Negative Binomial

void cdfnbn(int *which,double *p,double *q,double *s,double *xn,
            double *pr,double *ompr,int *status,double *bound)

  1     P and Q        [0,1];(0,1]       -----------        SUM to 1.0
  2     S              [0,infinity)      [0,1E100]
  3     XN             [0,infinity)      [0,1E100]
  4     PR and OMPR    [0,1];[0,1]       [0,1];[0,1]        SUM to 1.0

-------------------- Normal

void cdfnor(int *which,double *p,double *q,double *x,
            double *mean,double *sd,int *status,double *bound)

  1     P and Q        (0,1],(0,1]       -----------        SUM to 1.0
  2     X              (-inf.,inf.)      -----------
  3     MEAN           (-inf.,inf.)      -----------
  4     SD             (0,infinity)      -----------

-------------------- Poisson

void cdfpoi(int *which,double *p,double *q,double *s,
            double *xlam,int *status,double *bound)

  1     P and Q        [0,1],(0,1]       -----------        SUM to 1.0
  2     S              [0,infinity)      [0,1E100]
  3     XLAM           [0,infinity)      [0,1E100]

-------------------- Student's t

void cdft(int *which,double *p,double *q,double *t,double *df,
          int *status,double *bound)

  1     P and Q        (0,1],(0,1]       -----------        SUM to 1.0
  2     T              (-inf.,inf.)      [-1E100,1E100]
  3     DF             (0,infinity)      [1E-100,1E10]

-------------------- Noncentral t

void cdftnc(int *which,double *p,double *q,double *t,double *df,
            double *pnonc,int *status,double *bound)

  1     P and Q        [0,1-1E-16],              -----------
  2     T              [-infinity,infinity]      [-1E100,1E100]
  3     DF             (0,infinity)              [1E-100,1E10]
  4     PNONC          [-infinity,infinity)      [-1e4,1E4]
