#ifndef SAMPLING_H
#define SAMPLING_H
#define LENS 4096
#include "NGSNGS_cli.h"
#include "fasta_sampler.h"

struct Parsarg_for_Sampling_thread{
  fasta_sampler *reffasta;
  int threadno;
  int totalThreads; //<- this contains the total number of threads. This is shared among all threads

  int* FragLen;
  double* FragFreq;
  int No_Len_Val;
  double distparam1;
  double distparam2;
  int LengthType;

  char *NtQual_r1;
  char *NtQual_r2;
  ransampl_ws ***QualDist_r1;
  ransampl_ws ***QualDist_r2;

  double* MisMatch;
  int doMisMatchErr;
  int MisLength;

  double *NtErr_r1;
  double *NtErr_r2;

  int threadseed;
  size_t reads;
  size_t BufferLength;
  
  BGZF **bgzf_fp;

  samFile *SAMout;
  sam_hdr_t *SAMHeader;
  bam1_t **list_of_reads;
  int LengthData;
  int MaximumLength;

  int AddAdapt;
  const char* Adapter_1;
  const char* Adapter_2;

  int DoBriggs;
  int DoBriggsBiotin;
  float *BriggsParam;
  
  const char* SizeFile;
  const char* SizeFileFlag;
  int FixedSize;

  int maxreadlength;
  
  outputformat_e OutputFormat;
  seqtype_e SeqType;
  const char* QualFlag;

  int DoSeqErr;
  int NoAlign;
  char PolyNt;

  int rng_type;
};

void* Sampling_threads(void *arg);


#endif
