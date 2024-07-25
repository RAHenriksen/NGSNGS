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
#include <htslib/sam.h>

struct struct_for_amplicon_threads{
  BGZF* amplicon_in_fp;
  BGZF* amplicon_out_fp;  // Include the output file in the thread struct
  samFile* amplicon_in_sam;
  samFile* amplicon_out_sam; 
  bam_hdr_t *hdr;
  int filetype;
  int startLine;
  int endLine;
  int threadid;
  char* Briggs;
  char* Indel;
  int rng_type;
  long int Seed;
  
  int fixqual;
  int DoSeqErr;
  int DoSeqErrDist;
  ransampl_ws ***QualDistProfile;
  char *NtQual;
  double *NtErr;

  double* MisMatch;
  int MisLength;
  int doMisMatchErr;

  size_t moduloread;
  size_t totalLines;
  size_t linesPerThread;

  ampliconformat_e OutFormat;
};

void* ProcessFAFQ(void* args);

void* ProcessBAM(void* args);

void* AmpliconThreadInitialize(ampliconformat_e OutFormat,const char* Amplicon_in_fp,int filetype,const char* Amplicon_out_fp,const char* Subprofile,int threads,char *Briggs,char *Indel,int seed,int rng_type,int fixqual,const char* QualProfile,size_t moduloread,size_t totalLines,size_t linesPerThread);

#endif
