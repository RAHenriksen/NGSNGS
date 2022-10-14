/*
  part of NGSNGS
  program models briggs errors found in aDNA. Method assumes that seqlen>30 and that seq is 5->3
*/

#ifndef BRIGGS_H
#define BRIGGS_H
#include "mrand.h"
#include "Briggs.h"
/*
  orginal is a contiguous subsequence of our reference genome which we will call forward or + strand.
  We also assume that original is 5 to 3 prime.
  We also assume that original is 5' to 5' (oppositite strand).
  
            *********orginal************
  strand1   +++++++++++ MAKLE DRADSFADSF
  strand2
  In words this means that our original char* is the single strand representaiton of the double stranded fragment such that our original represents the 5 to 5.
  
  The output of the algorithm will be put in the results ** which is a char*[4] type. Each of these are assumed to be preallocated. The code in this function call does not allocate new memory and the caller are resposiable for cleaning up original and resutls.

 */


int SimBriggsModel2(char *original, int L, double nv, double lambda, double delta_s, double delta,mrand_t *mr,char **results);

#endif /* NGSNGSFUNC_H */
