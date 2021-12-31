#ifndef SAMPLING_H
#define SAMPLING_H

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,const char* Adapter_2,const char* OutputFormat,const char* SeqType);

#endif