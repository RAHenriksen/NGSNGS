/*
  part of NGSNGS
  program models briggs errors found in aDNA. Method assumes that seqlen>30 and that seq is 5->3
*/

#ifndef BRIGGS_H
#define BRIGGS_H
#include "mrand.h"

int Random_geometric_k(const double p,mrand_t *mr);

void SimBriggsModel(char seq[], int L, double nv, double lambda, double delta_s, double delta,mrand_t *mr,int rng_type);

#endif /* NGSNGSFUNC_H */
