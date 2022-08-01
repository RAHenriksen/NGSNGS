#ifndef THREADGENERATION_H
#define THREADGENERATION_H
#include "NGSNGS_cli.h"
#define LENS 4096

void* ThreadInitialization(faidx_t *seq_ref,int thread_no, int seed, size_t reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,outputformat_e OutputFormat,seqtype_e SeqType,float BriggsParam[4],const char* Briggs_Flag,
                        const char* Sizefile,int FixedSize,int SizeDistType,double val1,double val2,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,int ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* SubFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat, char* Variant_flag,const char *VarType,
                        char CommandArray[1024],const char* version,const char* HeaderIndiv,const char* NoAlign,size_t BufferLength);

#endif
