#ifndef AMPLICON_H
#define AMPLICON_H
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include "Briggs.h"
#include "mrand.h"
#include "fasta_sampler.h"
#include "NtSubModels.h"
#include "add_indels.h"
#include "HelpPage.h"
#include "Amplicon_cli.h"
#include <htslib/kstring.h>

struct struct_for_amplicon_threads{
  BGZF* amplicon_in_fp;
  int filetype;
  int startLine;
  int endLine;
  int threadid;
  char* Briggs;
  char* Indel;
  int rng_type;
  BGZF* amplicon_out_fp;  // Include the output file in the thread struct
  long int Seed;

  double* MisMatch;
  int MisLength;
  int doMisMatchErr;

  size_t moduloread;
  size_t totalLines;
  size_t linesPerThread;

  ampliconformat_e OutFormat;
};

void* AmpliconFAFQThreadInitialize(ampliconformat_e OutFormat,const char* Amplicon_in_fp,int filetype,const char* Amplicon_out_fp,const char* Subprofile,int threads,char *Briggs,char *Indel,int seed,int rng_type,size_t moduloread,size_t totalLines,size_t linesPerThread);

#endif
