#ifndef BRIGGS_H
#define BRIGGS_H
#include "mrand.h"

int Random_geometric_k(unsigned int  seed,const double p);

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed,mrand_t *mr);

#endif /* NGSNGSFUNC_H */