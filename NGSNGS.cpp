#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "Sampling.h"
#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "version.h"
#include "ThreadGeneration.h"
#include "HelpPage.h"
#include "NGSNGS.h"
#include "NGSNGS_cli.h"
#include "fasta_sampler.h"

#include <signal.h>
#define LENS 10000
#define MAXBINS 100

int VERBOSE = 1;
int SIG_COND = 1;

int really_kill =3;

void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
    fprintf(stderr,"\n\t-> If you really want ngsngs to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(1);
  VERBOSE=0;
  SIG_COND=0;
}

void catchkill(){
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  
}

int main(int argc,char **argv){
  // NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA - Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen 
  
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0))){
    fprintf(stderr,"\t-> ngsngs version %s: %s (htslib: %s) build(%s %s)\n",NGSNGS_RELEASE,NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    return 0;
  }
  else{
    catchkill();
    mypars = getpars(argc,argv);
    if(mypars==NULL)
      return 1;
    
    fprintf(stderr,"\n\t-> ngsngs version %s: %s (htslib: %s) build(%s %s)\n",NGSNGS_RELEASE,NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stderr,"\t-> My commmand: %s\n",mypars->CommandRun);

    clock_t t = clock();
    time_t t2 = time(NULL);

    outputformat_e OutputFormat = mypars->OutFormat;
    double readcov = mypars->coverage;

    // OS specific random number generator
    if (mypars->rng_type == -1){
      #if defined(__linux__) || defined(__unix__)
        mypars->rng_type = 0;
      #elif defined(__APPLE__) || defined(__MACH__)
        mypars->rng_type = 3;
        //when 0 it will have problems with drand48 reentrant, will default to erand48 (MacroRandType 3)
      #else
      #   error "Unknown compiler"
      #endif
    }
    
    if (mypars->Reference == NULL){ErrMsg(1.0);}
    if (mypars->OutName == NULL){ErrMsg(8.0);}
    
    int readcycle = 0;
    int nlines = 0;
    const char* SizeDist = mypars->LengthDist;
    double mean_length = 0;
    int SizeDistType=-1;double SizeDistval1 = 0; double SizeDistval2  = 0;
    
    // quality profile and cycle length relevant for specific formats (fastq and sequence alignment map format)
    if(OutputFormat==fqT|| OutputFormat== fqgzT ||OutputFormat==samT ||OutputFormat==bamT|| OutputFormat== cramT){
      if (mypars->CycleLength != 0){
        readcycle = mypars->CycleLength;
        if (mypars->QualProfile1 == NULL && mypars->FixedQual == 0){
          ErrMsg(11.0);
        }

        if(mypars->QualProfile1 != NULL){
          gzFile gz = Z_NULL;
          assert(((gz = gzopen(mypars->QualProfile1,"rb")))!=Z_NULL && "Check the structure of the provided nucleotide quality profile, see README on https://github.com/RAHenriksen/NGSNGS");
          gzclose(gz);
        }
      }
      else{
        if (mypars->QualProfile1 == NULL && mypars->FixedQual == 0){
          ErrMsg(11.0);
        }
        else if (mypars->QualProfile1 != NULL){
          gzFile gz = Z_NULL;
          assert(((gz = gzopen(mypars->QualProfile1,"rb")))!=Z_NULL && "Check the structure of the provided nucleotide quality profile, see README on https://github.com/RAHenriksen/NGSNGS");
          char buf[LENS];
          while(gzgets(gz,buf,LENS))
            nlines++;
          gzclose(gz);
          readcycle = (nlines-2)/5;
        }
      }
    }

    // Fragment length (-l, -lf, -ld)
    if (mypars->Length != 0){
      if (mypars->Length < 0){ErrMsg(3.0);}
      else{
        mean_length = mypars->Length;
        //fprintf(stderr,"Mean_length mypars-length %f\n",mean_length);
        SizeDistType=0;
        fprintf(stderr,"\t-> Mean fragment length equals the fixed length (-l) of %f nt\n",mean_length);
      } 
    }
    else if (mypars->LengthFile != NULL){
      FILE *fp = fopen(mypars->LengthFile, "r");
      int Len_prev; int Len_cur;
      float CummFreq_cur, CummFreq_prev;

      if (OutputFormat==faT|| OutputFormat== fagzT) {
        // Calculate mean length for all lines, i.e. the fragment length determines the number of reads, given a specific coverage
        if (fscanf(fp, "%d %f", &Len_prev, &CummFreq_prev) != 2) {
          // Handle error: unable to read expected data
          fclose(fp);
          return 1; // Or appropriate error code
        }
        mean_length += Len_prev * ((double)(CummFreq_prev));
        while (fscanf(fp, "%d %f", &Len_cur, &CummFreq_cur) == 2) {
          mean_length += Len_cur * ((double)(CummFreq_cur-CummFreq_prev));
          Len_prev = Len_cur;
          CummFreq_prev = CummFreq_cur;
        }
      }
      else{
        // Calculate mean length up to a certain readcycle limit, i.e. the read length determines the number of reads, given a specific coverage
        if (fscanf(fp, "%d %f", &Len_prev, &CummFreq_prev) != 2) {
          // Handle error: unable to read expected data
          fclose(fp);
          return 1; // Or appropriate error code
        }
        
        int Len_cur_tmp;float CummFreq_cur_tmp;
        if(mypars->FixedQual > 0 && mypars->CycleLength == 0){
          FILE *fp_tmp = fopen(mypars->LengthFile, "r");
          while (fscanf(fp_tmp, "%d %f", &Len_cur_tmp, &CummFreq_cur_tmp) == 2){continue;}
          fclose(fp_tmp);
          readcycle = (int) Len_cur_tmp;
        }
        
        //std::cout << " s 1 " <<mean_length << " z " << Len_prev << std::endl;
        mean_length += Len_prev * ((double)(CummFreq_prev));
        //std::cout << " s 2 " <<mean_length << std::endl;


        while (fscanf(fp, "%d %f", &Len_cur, &CummFreq_cur) == 2) {
          if (Len_prev < readcycle && Len_cur > readcycle) {
            // If the next line exceeds the limit, only consider the fraction of the current line up to the limit
            double fraction = (CummFreq_prev + (CummFreq_cur-CummFreq_prev) *
                              ((double)(readcycle-Len_prev)/(Len_cur-Len_prev)));
            mean_length += readcycle * fraction;
            break;
          }
          else if (Len_cur <= readcycle){
            // If the current line does not exceed the limit, add its full length to the mean
            mean_length += Len_cur * ((double)(CummFreq_cur-CummFreq_prev));
            Len_prev = Len_cur;
            CummFreq_prev = CummFreq_cur;
          }
        }
      }   
      fclose(fp);
      fprintf(stderr,"\t-> Mean fragment length of the provided length file (-lf) is %f nt\n",mean_length);
      if (mypars->Length <0){fprintf(stderr,"Fixed fragment length %d",mypars->Length);ErrMsg(5.0);}
      SizeDistType=1;
    }
    else if (SizeDist != NULL){
      char* Dist;
      char* DistParam = strdup(SizeDist);
      Dist = strtok(DistParam,",");
      SizeDistval1 = atof(strtok (NULL, ","));
      char* tmp = strtok(NULL, ",");
      if(tmp == NULL){SizeDistval2 = 0;}
      else{SizeDistval2 = atof(tmp);}
      
      if (strcasecmp(Dist,"Uni")==0){SizeDistType=2;mean_length=(double)(0.5*(SizeDistval1+SizeDistval2));}
      else if (strcasecmp(Dist,"Norm")==0){SizeDistType=3;mean_length= (double)SizeDistval1;}
      else if (strcasecmp(Dist,"LogNorm")==0){SizeDistType=4;mean_length= (double)exp((SizeDistval1+((SizeDistval2*SizeDistval2)/2)));}
      else if (strcasecmp(Dist,"Pois")==0){SizeDistType=5;mean_length= (double)SizeDistval1;}
      else if (strcasecmp(Dist,"Exp")==0){SizeDistType=6;mean_length= (double)(1/SizeDistval1);}
      else if (strcasecmp(Dist,"Gam")==0){SizeDistType=7;mean_length= (double)(SizeDistval1*SizeDistval2);}
      else if (mypars->Length >0){ErrMsg(5.0);}
      free((char *)Dist);
      fprintf(stderr,"\t-> Mean fragment length calculated from the provided length distribution (-ld) and parameters is %f nt\n",mean_length);
    }

    if (mypars->CycleLength == 0 && mypars->QualProfile1 == NULL){
      readcycle = (int) mean_length;
      fprintf(stderr,"\t-> The read cycle length is equal to the mean length from the length distribution (-ld) or fixed length (-l): %d\n",readcycle);
    }
    else{
      fprintf(stderr,"\t-> The read cycle length is either provided (-cl): %d or inferred from the quality profile dimension (-q1): %d\n",mypars->CycleLength,readcycle);
    }

    // read in reference genome to calculate the number of reads or coverage depending on genome length
    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(mypars->Reference);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);

    //first capture the cases where no reads or cov has been defined or both defined
    if(mypars->nreads == 0 && readcov == 0.0){
      fprintf(stderr,"must suply number of reads (-r) or desired coverage (-c)");exit(1);
    }
    if(mypars->nreads > 0 &&readcov > 0.0){
      fprintf(stderr,"must not suply number of reads (-r) and desired coverage (-c)");exit(1);
    }
    //fprintf(stderr,"NOW IM AFTER NREADS %zu \n",mypars->nreads);
    
    //now compute the number of reads required across all threads depending on genome size when providing a depth of coverage, NB. when -c > 1 the breadth of coverage will also converge towards 1 due to uniform sampling of positions 
    if (readcov > 0.0){
      size_t genome_size = 0;

      for (int i = 0; i < chr_total; i++){
        const char *chr_name = faidx_iseq(seq_ref,i);
        int chr_len = faidx_seq_len(seq_ref,chr_name);
        genome_size += chr_len;
      }
      if(OutputFormat==fqT|| OutputFormat== fqgzT ||OutputFormat==samT ||OutputFormat==bamT|| OutputFormat== cramT){
        if (mean_length >= readcycle){mean_length = (double)readcycle;}
      }
      if (mypars->seq_type == PE){
        mypars->nreads = ((readcov*genome_size)/mean_length)/2;
      }
      else{
        mypars->nreads = (readcov*genome_size)/mean_length;
      }
    }


    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in input file (-i): \'%s\': %d\n",mypars->Reference,chr_total);
    if(mypars->Glob_seed_binary == 1){
      fprintf(stderr,"\t-> Seed is provided (-s): %d\n",mypars->Glob_seed/1000);
    }
    else{
      fprintf(stderr,"\t-> Seed is set at current calendar time: %d\n",mypars->Glob_seed);
    }
    fprintf(stderr,"\t-> Number of sampling threads used (-t): %d and number of compression threads (-t2): %d\n",mypars->SamplThreads,mypars->CompressThreads);

    if (readcov > 0.0){fprintf(stderr,"\t-> Number of reads to be simulated %zu, based on the desired coverage (-c): %f\n",mypars->nreads,mypars->coverage);}
    else{fprintf(stderr,"\t-> Number of reads to be simulated (-r): %zu\n",mypars->nreads);}

    // additional nucleotides for sequencing reads 
    int AddAdapt = 0;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      AddAdapt = 1;

      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
      else{Polynt = "F";}
    }
    else{
      if (mypars->Poly != NULL){ErrMsg(14.0);}
      else{Polynt = "F";}
    }
    
    // read in quality profiles
    const char* QualStringFlag;
    if (mypars->QualProfile1 == NULL && mypars->FixedQual == 0){QualStringFlag = "false";}
    else{QualStringFlag = "true";}
    
    if (strcasecmp("true",QualStringFlag)==0){
      if(OutputFormat==fqT|| OutputFormat== fqgzT ||OutputFormat==samT ||OutputFormat==bamT|| OutputFormat== cramT){
        if (mypars->seq_type == PE && mypars->QualProfile2 == NULL && mypars->FixedQual == 0){
          //fprintf(stderr,"OUTPUTFORMAT 1");
          ErrMsg(11.0);
        }
      }
    }
    else{
      if(OutputFormat== fqT ||OutputFormat==fqgzT){
        //fprintf(stderr,"OUTPUTFORMAT 2");
        ErrMsg(11.0);
      }
    }

    //NB!
    /*if (strcasecmp("false",QualStringFlag)==0){
      if(strcasecmp("true",Adapt_flag==0 && mypars->Poly != NULL){WarMsg(2.0);}
      //if(mypars->ErrorFlag == NULL){WarMsg(3.0);}
    }
    if (strcasecmp("true",QualStringFlag)==0){
      if(OutputFormat== faT|| fagzT==OutputFormat)
	    WarMsg(4.0);
    }*/

    int qualstringoffset = 0;
    if(fqT==OutputFormat|| fqgzT==OutputFormat)
      qualstringoffset = 33;
    
    // PMD
    int DoNonBiotin = 0;
    int DoBiotin = 0;
    float PMDParam[4];
    if (mypars->PMD_NonBiotin != NULL || mypars->PMD_Biotin != NULL){
      char* BriggsParam;
      if (mypars->PMD_NonBiotin != NULL){
        BriggsParam = strdup(mypars->PMD_NonBiotin);
        DoNonBiotin = 1;
      }
      else{
        BriggsParam = strdup(mypars->PMD_Biotin);
        DoBiotin = 1;
      }
      
      PMDParam[0] = myatof(strtok(BriggsParam,"\", \t"));
      PMDParam[1] = myatof(strtok(NULL,"\", \t"));
      PMDParam[2] = myatof(strtok(NULL,"\", \t"));
      PMDParam[3] = myatof(strtok(NULL,"\", \t"));
      
      free(BriggsParam);
    }
    
    // Stochastic indels
    int DoIndel = 0;
    float IndelFuncParam[4];
    if (mypars->Indel != NULL){
      char* IndelInputParam = strdup(mypars->Indel);
      IndelFuncParam[0] = myatof(strtok(IndelInputParam,"\", \t"));
      IndelFuncParam[1] = myatof(strtok(NULL,"\", \t"));
      IndelFuncParam[2] = myatof(strtok(NULL,"\", \t"));
      IndelFuncParam[3] = myatof(strtok(NULL,"\", \t"));
      DoIndel = 1;
      free(IndelInputParam); 
    }

    int doMisMatchErr = 0;
    int MisMatchLength = 0;
    
    if (mypars->SubProfile != NULL)
      doMisMatchErr = 1;
    
    if (mypars->MisMatchMatrix_bdam != NULL)
      doMisMatchErr = 2;
    
    if(mypars->SubProfile != NULL && mypars->PMD_NonBiotin != NULL){
      ErrMsg(12.0);
    }
    else if(mypars->SubProfile != NULL && mypars->MisMatchMatrix_bdam != NULL){
      ErrMsg(12.0);
    }
    else if(mypars->MisMatchMatrix_bdam != NULL && mypars->PMD_NonBiotin != NULL){
      ErrMsg(12.0);
    }


    if(DoBiotin==1){
      fprintf(stderr,"\t-> Number of PCR duplicates for non-biotin simulated deamination model is %d\n",mypars->Duplicates);
    }
    
    /*
    int qsreadcyclelen = 0;
    if(OutputFormat==fqT|| OutputFormat== fqgzT ||OutputFormat==samT ||OutputFormat==bamT|| OutputFormat== cramT){      
      if (mypars->CycleLength > 0){
        qsreadcyclelen=1;
      }

      if (qsreadcyclelen == 0){
        fprintf(stderr,"\t-> When providing a fixed quality score (-qs) the read cycle length is not used as an upper limit for the simulated reads, unless -cl is specified\n");
      }
    }*/

    // Initialize the sampling threads
    ThreadInitialization(
      NGSNGS_VERSION,mypars->CommandRun,mypars->SamplThreads,mypars->CompressThreads,mypars->KstrBuf,
      mypars->Reference,mypars->Chromosomes,mypars->Glob_seed,mypars->rng_type,
      mypars->OutName,mypars->OutFormat,mypars->Align,mypars->seq_type,mypars->simmode,mypars->nreads, mypars->flankingregion,mypars->BedFile,mypars->MaskBed,
      mypars->LengthFile,mypars->Length,SizeDistType,SizeDistval1,SizeDistval2,mypars->LowerLimit,
      AddAdapt,mypars->Adapter1,mypars->Adapter2,Polynt,
      mypars->DoSeqErr,QualStringFlag,qualstringoffset,mypars->QualProfile1,mypars->QualProfile2,mypars->FixedQual,readcycle,mypars->CycleLength,
      doMisMatchErr,mypars->SubProfile,MisMatchLength,mypars->MisMatchMatrix_bdam,mypars->M3outname,
      PMDParam,DoNonBiotin,DoBiotin,mypars->Duplicates,
      mypars->mutationrate,mypars->referencevariations,mypars->generations,mypars->VariantDumpFile,
      mypars->vcffile,mypars->HeaderIndiv,mypars->NameIndiv,mypars->VCFDumpFile,mypars->CaptureVCF,mypars->linkage,
      IndelFuncParam,DoIndel,mypars->IndelDumpFile);

    /*
    NGSNGS overall - @param const char* version,char CommandArray[LENS],int thread_no,int threadwriteno,size_t BufferLength
    Reference specific - @param const char* refSseq,const char* Specific_Chr,int seed,int RandMacro,
    simulation and output specific - @param const char* OutputName,outputformat_e OutputFormat,int Align,seqtype_e SeqType,int simmode,size_t reads,size_t flankingregion, const char* BedFile, int MaskBed,
    fragment length specific - @param const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int Lowerlimit,
    Additional nucleotide post simulation specific - @param int AddAdapt,const char* Adapter_1,const char* Adapter_2,const char* Polynt,
    Sequencing error (fastq,sam,bam,cram) specific - @param int DoSeqErr,const char* QualStringFlag,int qualstringoffset,const char* QualProfile1,const char* QualProfile2,int FixedQual,int readcycle,int readcycle_fix,
    Nucleotide misincorporation specific - @param int doMisMatchErr,const char* SubProfile,int MisLength,const char* MisMatchMatrix,const char* M3outname,
    PMD specific - @param float BriggsParam[4],int DoNonBiotin,int DoBiotin,int Duplicates,
    Reference specific stochastic variation - @param double mutationrate, size_t referencevariations, int generations,
    Allele specific variations - @param const char *VariantFile,int HeaderIndivIdx,const char* NameIndiv,const char* VCFfileDump,int CaptureVCF,int linkage,
    sequencing read specific stochastic indels variations - @param float IndelFuncParam[4],int DoIndel,const char* IndelDumpFile
    */

    fai_destroy(seq_ref);
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  argStruct_destroy(mypars);
}

