#ifndef FASTA_SAMPLER_H
#define FASTA_SAMPLER_H
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include <htslib/faidx.h>//for faidx
#include "RandSampling.h"
#include "mrand.h"

typedef struct{
  faidx_t *fai;
  int nref;
  char **seqs;
  char **seqs_names;
  int *seqs_l;
  size_t seq_l_total;
  ransampl_ws *ws;
}fasta_sampler;


fasta_sampler *fasta_sampler_alloc(const char *,const char *);
void fasta_sampler_destroy(fasta_sampler *fs);
char* sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength);
#endif
