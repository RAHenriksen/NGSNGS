#ifndef ADD_INDEL_H
#define ADD_INDEL_H
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include "mrand.h"
#include "fasta_sampler.h"
#include <htslib/kstring.h>
#include "RandSampling.h"

//pars is a general placeholder used for representing general parameters. Currently it is , insertion rate (prob), deletion rate (prob), insertion length(expontial lambda),deletion length(expontial lambda),
void add_indel(mrand_t *mr,char *fragmentn,int readlength,double *pars,char *INDEL_INFO,int *ops);

//void add_indel_kstring(mrand_t *mr,kstring_t* seq,double *pars,char *INDEL_INFO,int* ops); As of 27-08-2024 it is to be added 

void add_indel_amplicon_fa(mrand_t *mr,kstring_t* seq,double *pars,int* ops);

void add_indel_amplicon_fqbam(mrand_t *mr,kstring_t* seq,kstring_t* qual,double *pars,int* ops,int ErrProbTypeOffset);

#endif
