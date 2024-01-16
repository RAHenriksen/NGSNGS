#ifndef SAMPLING_H
#define SAMPLING_H
#define LENS 10000
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
  int FixedQual_r1r2;
  int qsreadcycle;
  
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
  int Duplicates;
  
  const char* SizeFile;
  const char* SizeFileFlag;
  int FixedSize;

  int lowerlimit;
  int maxreadlength;
  
  seqtype_e SeqType;
  const char* QualFlag;

  int DoSeqErr;
  int Align;
  char PolyNt;

  int rng_type;
  float *IndelFuncParam;
  int DoIndel;
  const char* IndelDumpFile;
};

void* Sampling_threads(void *arg);


#endif
