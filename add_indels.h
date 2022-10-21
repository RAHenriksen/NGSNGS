#ifndef ADD_INDEL_H
#define ADD_INDEL_H
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include "mrand.h"
#include "fasta_sampler.h"

int Random_geometric_k(const double p,mrand_t *mr);
//pars is a general placeholder used for representing general parameters. Currently it is , insertion rate (prob), length(expontial lambda),deletion rate (prob), length(expontial lambda),
int add_indel(mrand_t *mr,char *fragmentn,int readlength,double *pars);
void add_indel_fs(fasta_sampler *fs,mrand_t *mr,double *pars);
#endif
