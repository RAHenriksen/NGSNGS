#ifndef THREADGENERATION_H
#define THREADGENERATION_H
#include "NGSNGS_cli.h"
#define LENS 10000

  /*
  NGSNGS overall - @param const char* version,char CommandArray[LENS],int thread_no,int threadwriteno,size_t BufferLength
  Reference specific - @param const char* refSseq,const char* Specific_Chr,int seed,int RandMacro,
  simulation and output specific - @param const char* OutputName,outputformat_e OutputFormat,int Align,seqtype_e SeqType,int simmode,size_t reads,size_t flankingregion, const char* BedFile, int MaskBed,
  fragment length specific - @param const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int Lowerlimit,
  Additional nucleotide post simulation specific - @param int AddAdapt,const char* Adapter_1,const char* Adapter_2,const char* Polynt,
  Sequencing error (fastq,sam,bam,cram) specific - @param int DoSeqErr,const char* QualStringFlag,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int FixedQual,int readcycle,int readcycle_fix,
  Nucleotide misincorporation specific - @param int doMisMatchErr,const char* SubProfile,int MisLength,const char* MisMatchMatrix,const char* M3outname,
  PMD specific - @param float BriggsParam[4],int DoBriggs,int DoBriggsBiotin,int Duplicates,
  Reference specific stochastic variation - @param double mutationrate, size_t referencevariations, int generations,
  Allele specific variations - @param const char *VariantFile,int HeaderIndivIdx,const char* NameIndiv,const char* VCFfileDump,int CaptureVCF,int linkage,
  sequencing read specific stochastic indels variations - @param float IndelFuncParam[4],int DoIndel,const char* IndelDumpFile
  */
  
void* ThreadInitialization(const char* version,char CommandArray[LENS],int thread_no,int threadwriteno,size_t BufferLength,
                        const char* refSseq,const char* Specific_Chr,int seed,int RandMacro,
                        const char* OutputName,outputformat_e OutputFormat,int Align,seqtype_e SeqType,int simmode,size_t reads,size_t flankingregion, const char* BedFile, int MaskBed,
                        const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int Lowerlimit,
                        int AddAdapt,const char* Adapter_1,const char* Adapter_2,const char* Polynt,
                        int DoSeqErr,const char* QualStringFlag,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int FixedQual,int readcycle,int readcycle_fix,
                        int doMisMatchErr,const char* SubProfile,int MisLength,const char* MisMatchMatrix,const char* M3outname,
                        float BriggsParam[4],int DoBriggs,int DoBriggsBiotin,int Duplicates,
                        double mutationrate, size_t referencevariations, int generations,char* VariationfileDump,
                        const char *VariantFile,int HeaderIndivIdx,const char* NameIndiv,const char* VCFfileDump,int CaptureVCF,int linkage,
                        float IndelFuncParam[4],int DoIndel,const char* IndelDumpFile);

#endif
