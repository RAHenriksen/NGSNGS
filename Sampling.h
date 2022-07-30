#ifndef SAMPLING_H
#define SAMPLING_H
#define LENS 4096

struct Parsarg_for_Sampling_thread{
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  char *genome;
  int *chr_idx_array;
  int chr_no;
  int threadno;
  int totalThreads; //<- this contains the total number of threads. This is shared among all threads
  size_t *size_cumm;
  char **names;

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
  const char* SubFlag;
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

  const char* Adapter_flag;
  const char* Adapter_1;
  const char* Adapter_2;

  const char* Briggs_flag;
  float *BriggsParam;
  
  const char* SizeFile;
  const char* SizeFileFlag;
  int FixedSize;

  int readcycle;
  
  const char* OutputFormat;
  const char* SeqType;
  const char* QualFlag;

  char ErrorFlag;
  char* Variant_flag;
  char NoAlign;
  char PolyNt;

  int rng_type;
};

void* Sampling_threads(void *arg);


#endif
