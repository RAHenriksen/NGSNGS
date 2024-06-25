/*
  part of NGSNGS
  program models briggs errors found in aDNA. Method assumes that seqlen>30 and that seq is 5->3
*/

#ifndef BRIGGS2_H
#define BRIGGS2_H
#include "mrand.h"
#include "Briggs2.h"
#include <htslib/kstring.h>

/*
  orginal is a contiguous subsequence of our reference genome which we will call forward or + strand.
  We also assume that original is 5 to 3 prime.
  We also assume that original is 5' to 5' (oppositite strand).
  
  In words this means that our original char* is the single strand representaiton of the double stranded fragment such that our original represents the 5 to 5.
  
  The output of the algorithm will be put in the results ** which is a char*[4] type. Each of these are assumed to be preallocated. The code in this function call does not allocate new memory and the caller are resposiable for cleaning up original and resutls.

 */


int SimBriggsModel2(char *original, int L, double nv, double lambda, double delta_s, double delta,mrand_t *mr,char **results, int strandR1,int& C_to_T_counter,int& G_to_A_counter,int& C_total,int& G_total);

int SimBriggsModel2_k(kstring_t* ori, int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr, kstring_t* res[], int strandR1,
                    int& C_to_T_counter, int& G_to_A_counter,
                    int& C_total, int& G_total);

#endif /* NGSNGSFUNC_H */
