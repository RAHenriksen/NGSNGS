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

#include "NGSNGS_func.h"
#include "Sampling.h"
#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "version.h"
#include "ThreadGeneration.h"
#include "HelpPage.h"
#include "NGSNGS.h"

#include <signal.h>
#define LENS 4096
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
    exit(0);
  VERBOSE=0;
  SIG_COND=0;
  //  pthread_mutex_unlock(&mUpPile_mutex);
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
  //fprintf(stderr,"PRINT TYPE %d\n",MacroRandType);
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    catchkill();
    mypars = getpars(argc,argv);
    if(mypars==NULL)
      return 1;
    
    char* Command = mypars->CommandRun;
    fprintf(stderr,"\n\t-> ngsngs version: %s (htslib: %s) build(%s %s)\n",NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stderr,"\t-> Mycommmand: %s\n",Command);

    //fprintf(stderr,"\t-> Command 2 : %s and version %s \n",CommandArray,version);
    clock_t t = clock();
    time_t t2 = time(NULL);
    int Glob_seed = mypars->Glob_seed; 

    const char *fastafile = mypars->Reference;
    const char* OutputFormat = mypars->OutFormat;
    const char* filename = mypars->OutName; //"chr22_out";
    const char* Seq_Type = mypars->Seq;
    //    size No_reads = mypars->nreads;
    double readcov = mypars->coverage;
    int MacroRandType = mypars->rand_val; //extern int

    if (MacroRandType == -1){
      #if defined(__linux__) || defined(__unix__) // all unices not caught above
      // Unix
        MacroRandType = 0;
      #elif defined(__APPLE__) || defined(__MACH__)
        MacroRandType = 3;
        //when 0 it will have problems with drand48 reentrant, will default to erand48 (MacroRandType 3)
      #else
      #   error "Unknown compiler"
      #endif
    }
    //fprintf(stderr,"RANDOM VALUE %d \n",MacroRandType);
    
    
    if (fastafile == NULL){ErrMsg(1.0);}
    if (Seq_Type == NULL){ErrMsg(6.0);}
    if (OutputFormat == NULL){ErrMsg(7.0);}
    if (filename == NULL){ErrMsg(8.0);}
    
    //fprintf(stderr,"\t-> Command: %s \n",Command);
    int FixedSize = mypars->Length;
    const char* Sizefile = mypars->LengthFile;
    const char* SizeDist = mypars->LengthDist;
    double meanlength = 0;
    int SizeDistType=-1;double val1 = 0; double val2  = 0;

    if (FixedSize != 0){
      if (FixedSize < 0){ErrMsg(3.0);}
      else{
        meanlength = FixedSize;
        SizeDistType=0;
      } 
    }
    if (Sizefile != NULL){
      //fprintf(stderr,"SIZE FILE ARG\n");
      double sum,n;
      sum=n=0;

      char buf[LENS];
      gzFile gz = Z_NULL;
      gz = gzopen(Sizefile,"r");
      assert(gz!=Z_NULL);
      while(gzgets(gz,buf,LENS)){
        double Length_tmp = atof(strtok(buf,"\n\t "));
        double Frequency_tmp = atof(strtok(NULL,"\n\t "));
        sum += Length_tmp*Frequency_tmp;
        n = n+1;
      }
      gzclose(gz);
      
      meanlength = sum/n;
      if (FixedSize <0){fprintf(stderr,"FIXED SIZE %d",FixedSize);ErrMsg(5.0);}
      SizeDistType=1;
    }
    if (SizeDist != NULL){
      char* Dist;
      char* DistParam = strdup(SizeDist);
      Dist = strtok(DistParam,",");
      val1 = atof(strtok (NULL, ","));
      char* tmp = strtok(NULL, ",");
      if(tmp == NULL){val2 = 0;}
      else{val2 = atof(tmp);}
      
      if (strcasecmp(Dist,"Uni")==0){SizeDistType=2;meanlength=(0.5*(val1+val2));}
      if (strcasecmp(Dist,"Norm")==0){SizeDistType=3;meanlength= val1;}
      if (strcasecmp(Dist,"LogNorm")==0){SizeDistType=4;meanlength= exp((val1+((val2*val2)/2)));}
      if (strcasecmp(Dist,"Pois")==0){SizeDistType=5;meanlength= val1;}
      if (strcasecmp(Dist,"Exp")==0){SizeDistType=6;meanlength= 1/val1;}
      if (strcasecmp(Dist,"Gam")==0){SizeDistType=7;meanlength= (val1/val2);}
      if (FixedSize >0){ErrMsg(5.0);}
      free((char *)Dist);
    }

    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int SamplThreads = mypars->SamplThreads;
    int CompressThreads = mypars->CompressThreads;

    //first capture the awkward cases where no reads or cov has been defined or both defined
    if(mypars->nreads == 0 && readcov == 0.0){
      fprintf(stderr,"must suply nreads or cov");exit(0);
    }
    if(mypars->nreads > 0 &&readcov > 0.0){
      fprintf(stderr,"must not suply nreads and cov");exit(0);
    }
    //now compute the number of reads required acroos all threads

    if (readcov > 0.0){
      size_t genome_size = 0;

      for (int i = 0; i < chr_total; i++){
        const char *chr_name = faidx_iseq(seq_ref,i);
        int chr_len = faidx_seq_len(seq_ref,chr_name);
        genome_size += chr_len;
      }
      mypars->nreads =  (readcov*genome_size)/meanlength;
      fprintf(stderr,"\t-> Number of simulated reads: %zu or coverage: %f\n",mypars->nreads,mypars->coverage);
    }
  
    size_t nreads_per_thread = mypars->nreads/SamplThreads;
    
    size_t BufferLength = mypars->KstrBuf;

    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",Glob_seed);
    fprintf(stderr,"\t-> Number of sampling threads used (-t): %d and number of compression threads (-t2): %d\n",SamplThreads,CompressThreads);
    fprintf(stderr,"\t-> Number of simulated reads: %zu or coverage: %f\n",mypars->nreads,mypars->coverage);

    const char* Adapt_flag;
    const char* Adapter_1 = NULL;
    const char* Adapter_2 = NULL;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      //fprintf(stderr,"\t-> ARGPARSE ADAPTER + POLY\n");
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      Adapter_2 = mypars->Adapter2;

      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
      else{Polynt = "F";}
      
    }
    else{
      //fprintf(stderr,"\t-> ARGPARSE ADAPT FLAG+ POLY\n");
      Adapt_flag = "false";
      if (mypars->Poly != NULL){ErrMsg(14.0);exit(0);}
      else{Polynt = "F";}
    }
    // QUALITY PROFILES
    const char* QualProfile1; const char* QualProfile2;
    QualProfile1 = mypars->QualProfile1; QualProfile2 = mypars->QualProfile2;

    const char* QualStringFlag;
    if (QualProfile1 == NULL){QualStringFlag = "false";}
    else{QualStringFlag = "true";}
    //fprintf(stderr,"qualstring test %s",QualStringFlag);
    if (strcasecmp("true",QualStringFlag)==0){
      if(OutputFormat && (strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0 || strcasecmp("sam",OutputFormat)==0 || strcasecmp("bam",OutputFormat)==0 || strcasecmp("cram",OutputFormat)==0)){
        if (strcasecmp("PE",Seq_Type)==0 && QualProfile2 == NULL){
          ErrMsg(11.0);
          //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
          exit(0);
        }
      }
    }
    else
    {
      if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){
        ErrMsg(11.0);
        //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
    }
    //fprintf(stderr,"\t-> ADAPTER FLAG IS:%s\n",Adapt_flag);
    //fprintf(stderr,"\t-> QUAL STRING FLAG IS:%s\n",QualStringFlag);
    const char* ErrorFlag;
    if (mypars->ErrorFlag != NULL){
      ErrorFlag = mypars->ErrorFlag;
    }
    else{ErrorFlag = "T";}
    
    const char* NoAlign;
    if (mypars->NoAlign != NULL){
      NoAlign = mypars->NoAlign;
    }
    else{NoAlign = "F";}

    if (strcasecmp("false",QualStringFlag)==0){
      if(strcasecmp("true",Adapt_flag)==0 && mypars->Poly != NULL){WarMsg(2.0);}
      if(mypars->ErrorFlag != NULL){WarMsg(3.0);}
    }
    if (strcasecmp("true",QualStringFlag)==0){
      if(strcasecmp("fa",OutputFormat)==0 || strcasecmp("fa.gz",OutputFormat)==0){WarMsg(4.0);}
    }

    int qualstringoffset = 0;
    if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){qualstringoffset = 33;}
    
    const char* Briggs_Flag;
    float Param[4];
    if (mypars->Briggs != NULL){
      char* BriggsParam = strdup(mypars->Briggs);
      Param[0] = myatof(strtok(BriggsParam,"\", \t"));
      Param[1] = myatof(strtok(NULL,"\", \t"));
      Param[2] = myatof(strtok(NULL,"\", \t"));
      Param[3] = myatof(strtok(NULL,"\", \t"));
      Briggs_Flag = "True";
      free(BriggsParam); // Again using strdup
    }
    else{Briggs_Flag = "False";}
    
    const char* SubProfile; const char* SubFlag;
    
    SubProfile = mypars->SubProfile;
    if (SubProfile == NULL){SubFlag = "false";}
    else{SubFlag = "true";}
    //fprintf(stderr,"SUB FLAG IS %s\n",SubFlag);
    if(SubProfile != NULL && mypars->Briggs != NULL){
      ErrMsg(12.0);
      exit(0);
    }

    const char* Specific_Chr[1024] = {};
    if (mypars->Chromosomes != NULL){
      fprintf(stderr,"PARTIAL chromosomes %s\n",mypars->Chromosomes);
      int chr_idx_partial = 0;
      Specific_Chr[chr_idx_partial++] = strtok(strdup(mypars->Chromosomes),"\", \t");
      char *chrtok = NULL;
      while(((chrtok=strtok(NULL,"\", \t")))){
	      Specific_Chr[chr_idx_partial++] = strdup(chrtok);
	      assert(chr_idx_partial<MAXBINS);
      }
      fprintf(stderr,"AFTER WHILE and chr_idx_partial %d\n",chr_idx_partial);
      Specific_Chr[chr_idx_partial++] = "\0";
    }
    
    char* Variant_flag =NULL;
    const char* VCFformat = mypars->Variant;
    const char* VarType =mypars->Variant_type;
    const char* HeaderIndiv = mypars->HeaderIndiv;

    if (VCFformat != NULL){
      Variant_flag = strdup("bcf");
      if(VarType == NULL){
        VarType = strdup("snp");
      }
      else if(VarType != NULL){
        VarType = mypars->Variant_type;
      }
      fprintf(stderr,"VARIANT TYPE %s\n",VarType);
    }

    //if(Specific_Chr[0]=='\0'){fprintf(stderr,"HURRA");}
    int DeamLength = 0;
    //const char* HeaderIndiv = "HG00096";
    fprintf(stderr,"LENGTH TYPE %d\n",SizeDistType);
    ThreadInitialization(seq_ref,SamplThreads,Glob_seed,nreads_per_thread,filename,
                      Adapt_flag,Adapter_1,Adapter_2,OutputFormat,Seq_Type,
                      Param,Briggs_Flag,Sizefile,FixedSize,SizeDistType,val1,val2,
                      qualstringoffset,QualProfile1,QualProfile2,CompressThreads,QualStringFlag,Polynt,
                      ErrorFlag,Specific_Chr,fastafile,SubFlag,SubProfile,DeamLength,MacroRandType,
                      VCFformat,Variant_flag,VarType,Command,NGSNGS_VERSION,HeaderIndiv,NoAlign,BufferLength);
    fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  // MEMORY DEALLOCATION OF STRDUP FROM INPUT PARAMETERS
  // REQUIRED DEALLOCATIONS
  free((char *)mypars->Reference); //-i
  free((char *)mypars->Seq); //
  
  free((char *)mypars->OutFormat);
  free((char *)mypars->OutName);
  free((char *)mypars->LengthFile);
  free((char *)mypars->LengthDist);
  
  free((char*)mypars->CommandRun);

  // OPTIONAL DEALLOCATIONS
  free((char *)mypars->Variant);
  free((char *)mypars->Variant_type);
  free((char *)mypars->Variant_type);
  free((char *)mypars->HeaderIndiv);
  free((char *)mypars->Adapter1);
  free((char *)mypars->Adapter2);
  free((char *)mypars->QualProfile1);
  free((char *)mypars->QualProfile2);
  free((char *)mypars->SubProfile);
  free((char *)mypars->Briggs);
  free((char *)mypars->Poly);
  free((char *)mypars->HeaderIndiv);
  delete mypars;
}