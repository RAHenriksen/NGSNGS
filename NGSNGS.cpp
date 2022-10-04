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
    
    fprintf(stderr,"\n\t-> ngsngs version: %s (htslib: %s) build(%s %s)\n",NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stderr,"\t-> Mycommmand: %s\n",mypars->CommandRun);

    clock_t t = clock();
    time_t t2 = time(NULL);

    //const char* OutputFormat = mypars->OutFormat;
    double readcov = mypars->coverage;

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
    //fprintf(stderr,"RANDOM VALUE %d \n",MacroRandType);
    
    if (mypars->Reference == NULL){ErrMsg(1.0);}
    if (mypars->OutName == NULL){ErrMsg(8.0);}
    
    const char* SizeDist = mypars->LengthDist;
    double meanlength = 0;
    int SizeDistType=-1;double val1 = 0; double val2  = 0;

    if (mypars->Length != 0){
      if (mypars->Length < 0){ErrMsg(3.0);}
      else{
        meanlength = mypars->Length;
        SizeDistType=0;
      } 
    }
    if (mypars->LengthFile != NULL){
      //fprintf(stderr,"SIZE FILE ARG\n");
      double sum,n;
      sum=n=0;

      char buf[LENS];
      gzFile gz = Z_NULL;
      gz = gzopen(mypars->LengthFile,"r");
      assert(gz!=Z_NULL);
      while(gzgets(gz,buf,LENS)){
        double Length_tmp = atof(strtok(buf,"\n\t "));
        double Frequency_tmp = atof(strtok(NULL,"\n\t "));
        sum += Length_tmp*Frequency_tmp;
        n = n+1;
      }
      gzclose(gz);
      
      meanlength = sum/n;
      if (mypars->Length <0){fprintf(stderr,"FIXED SIZE %d",mypars->Length);ErrMsg(5.0);}
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
      if (mypars->Length >0){ErrMsg(5.0);}
      free((char *)Dist);
    }

    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(mypars->Reference);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);

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
    }
  
    //size_t nreads_per_thread = mypars->nreads/mypars->SamplThreads;
    
    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",mypars->Reference,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",mypars->Glob_seed);
    fprintf(stderr,"\t-> Number of sampling threads used (-t): %d and number of compression threads (-t2): %d\n",mypars->SamplThreads,mypars->CompressThreads);
    fprintf(stderr,"\t-> Number of simulated reads: %zu or coverage: %f\n",mypars->nreads,mypars->coverage);

    int AddAdapt = 0;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      AddAdapt = 1;

      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
      else{Polynt = "F";}
    }
    else{
      //fprintf(stderr,"\t-> ARGPARSE ADAPT FLAG+ POLY\n");
      if (mypars->Poly != NULL){ErrMsg(14.0);exit(0);}
      else{Polynt = "F";}
    }
    // QUALITY PROFILES

    const char* QualStringFlag;
    if (mypars->QualProfile1 == NULL){QualStringFlag = "false";}
    else{QualStringFlag = "true";}
    //fprintf(stderr,"qualstring test %s",QualStringFlag);
    outputformat_e OutputFormat = mypars->OutFormat;
    if (strcasecmp("true",QualStringFlag)==0){
      if(OutputFormat==fqT|| OutputFormat== fqgzT ||OutputFormat==samT ||OutputFormat==bamT|| OutputFormat== cramT){
        if (mypars->seq_type == PE && mypars->QualProfile2 == NULL){
          ErrMsg(11.0);
          //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
          exit(0);
        }
      }
    }
    else
    {
      if(OutputFormat== fqT ||OutputFormat==fqgzT){
        ErrMsg(11.0);
        //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
    }
    //fprintf(stderr,"\t-> ADAPTER FLAG IS:%s\n",Adapt_flag);
    //fprintf(stderr,"\t-> QUAL STRING FLAG IS:%s\n",QualStringFlag);
    

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
    
    int DoBriggs = 0;
    int DoBriggsBiotin = 0;
    float Param[4];
    if (mypars->Briggs != NULL || mypars->BriggsBiotin != NULL){
      char* BriggsParam;
      if (mypars->Briggs != NULL){
        BriggsParam = strdup(mypars->Briggs);
        DoBriggs = 1;
      }
      else{
        BriggsParam = strdup(mypars->BriggsBiotin);
        DoBriggsBiotin = 1;
      }
      
      Param[0] = myatof(strtok(BriggsParam,"\", \t"));
      Param[1] = myatof(strtok(NULL,"\", \t"));
      Param[2] = myatof(strtok(NULL,"\", \t"));
      Param[3] = myatof(strtok(NULL,"\", \t"));
      
      free(BriggsParam); // Again using strdup
    }
    
    int doMisMatchErr = 0;
    
    if (mypars->SubProfile != NULL)
      doMisMatchErr = 1;

    if(mypars->SubProfile != NULL && mypars->Briggs != NULL){
      ErrMsg(12.0);
      exit(0);
    }

    int DeamLength = 0;
    //mypars->nreads/mypars->SamplThreads
    ThreadInitialization(mypars->Reference,mypars->SamplThreads,mypars->Glob_seed,mypars->nreads,mypars->OutName,
                      AddAdapt,mypars->Adapter1,mypars->Adapter2,mypars->OutFormat,mypars->seq_type,
                      Param,DoBriggs,DoBriggsBiotin,mypars->LengthFile,mypars->Length,SizeDistType,val1,val2,
                      qualstringoffset,mypars->QualProfile1,mypars->QualProfile2,mypars->CompressThreads,QualStringFlag,Polynt,
                      mypars->DoSeqErr,mypars->Chromosomes,doMisMatchErr,mypars->SubProfile,DeamLength,mypars->rng_type,
                      mypars->vcffile,mypars->CommandRun,NGSNGS_VERSION,mypars->HeaderIndiv,mypars->NoAlign,mypars->KstrBuf);
    fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  argStruct_destroy(mypars);
}
