#ifndef SAMPLING_H
#define SAMPLING_H
#define LENS 4096

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, size_t reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_Flag,
                        const char* Sizefile,int FixedSize,int SizeDistType,int val1,int val2,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,const char* ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* SubFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat, char* Variant_flag,const char *VarType,
                        char CommandArray[1024],const char* version,const char* HeaderIndiv,const char* NoAlign,size_t BufferLength);

#endif
