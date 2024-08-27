/*
  part of NGSNGS
  program models briggs errors found in aDNA. Method assumes that seqlen>30 and that seq is 5->3
*/

#ifndef BRIGGS_H
#define BRIGGS_H
#include "mrand.h"
#include <htslib/kstring.h>

int Biotin_ds_454Roche(char seq[], int L, double nv, double lambda, double delta_s, double delta,mrand_t *mr,int strand,int& C_to_T_counter,int& G_to_A_counter,int& C_to_T_counter_rev,int& G_to_A_counter_rev);

int Biotin_ds_454Roche_kstring(kstring_t* seq, double nv, double lambda, double delta_s, double delta, mrand_t* mr, int strand, int& C_to_T_counter, int& G_to_A_counter, int& C_to_T_counter_rev, int& G_to_A_counter_rev);
// the kstring equivalent is as of 27-08-2024 not incorporated into the actual NGSNGS functionality but is meant to be for long-read sequence simulation

int PMD_Amplicon(kstring_t* seq, double nv, double lambda, double delta_s, double delta, mrand_t* mr);

#endif /* NGSNGSFUNC_H */
