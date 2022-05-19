#ifndef SAMPLING_H
#define SAMPLING_H
#define LENS 4096

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_Flag,
                        const char* Sizefile,int FixedSize,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,const char* ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* SubFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat,const char* Variant_flag,const char *VarType);

#endif