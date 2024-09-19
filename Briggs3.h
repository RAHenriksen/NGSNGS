/*
  part of NGSNGS
  program models briggs errors found in aDNA. Method assumes that seqlen>30 and that seq is 5->3
*/

#ifndef BRIGGS3_H
#define BRIGGS3_H
#include "mrand.h"
#include "Briggs3.h"
#include <htslib/kstring.h>

/*
  orginal is a contiguous subsequence of our reference genome which we will call forward or + strand.
  We also assume that original is 5 to 3 prime.
  We also assume that original is 5' to 5' (oppositite strand).
  
  In words this means that our original char* is the single strand representaiton of the double stranded fragment such that our original represents the 5 to 5.
  
  The output of the algorithm will be put in the results ** which is a char*[4] type. Each of these are assumed to be preallocated. The code in this function call does not allocate new memory and the caller are resposiable for cleaning up original and resutls.

 */



int SimBriggsModel3(char *ori, int L, double nv, double lambda, double delta_s, double delta_d,double end3primeoverhang, mrand_t *mr,char **&res,int* &frag_type,int strandR1,
                    int& C_to_T_counter,int& G_to_A_counter,int& C_total,int& G_total,
                    int &total_frag,int& fwd_fragment_count,int& rev_fragment_count);

#endif /* BRIGGS3_H */
