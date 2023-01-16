#include "../mrand.h"
#include "../fasta_sampler.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <iostream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

typedef struct{
  size_t RandRefLen;
  int seed;
  int ChrNo;
  int DoNnt;
  const char* RandRef_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Generate artificial reference genome uniformly sampled from the four nucleotides\n");
  fprintf(fp,"Usage\n./RandRef -l <Genome Length> -mp <Modulus position> -o <Prefix of output name>\n");
  fprintf(fp,"\nExample\n\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-l   | --length: \t\t Reference genome .fa format\n");
  fprintf(fp,"-s   | --seed: \t\t\t seed for random generators\n");
  fprintf(fp,"-p   | --ploidy: \t\t Reference genome .fa format\n");
  fprintf(fp,"-n   | --Nt: \t\t Add refererence nucleotide 'N-nucleotide'\n");
  fprintf(fp,"-pos | --position: \t\t\t Internal txt file with chromomsomal coordinate for the incorporated variations\n");
  fprintf(fp,"-o   | --output: \t\t Altered reference genomes saved as .fa\n");
  exit(1);
  return 0;
}
 

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->RandRefLen = 0;
  mypars->seed = 0;
  mypars->ChrNo = 1;
  mypars->DoNnt = 0;
  mypars->RandRef_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      mypars->RandRefLen = atol(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->seed = atoi(*(++argv));
    }
    else if(strcasecmp("-n",*argv)==0 || strcasecmp("--chrNo",*argv)==0){
      mypars->ChrNo = atoi(*(++argv));
    }
    else if(strcasecmp("-nn",*argv)==0 || strcasecmp("--NoNt",*argv)==0){
      mypars->DoNnt = 1;
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->RandRef_out = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}


int Reference_Variant(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    else{
      mypars = getpars(argc,argv);
      size_t RandRefLen = mypars->RandRefLen;
      int seed = mypars->seed;
      int ChrNo = mypars->ChrNo;
      int DoNnt = mypars->DoNnt;
      const char* RandRef_out = mypars->RandRef_out;
      
      mrand_t *mr = mrand_alloc(0,seed);
      long rand_val;
      
      const char *bases; 
      int ntval;

      if(DoNnt == 0){
        bases = "ACGTN";
        ntval = 5;
      }
      else{
        bases = "ACGT";
        ntval = 4;
      }    
      
      BGZF **bgzf_fp = (BGZF **) calloc(1,sizeof(BGZF *));
      int mt_cores = 1;
      int bgzf_buf = 256;
      const char* mode = "wu";

      bgzf_fp[0] = bgzf_open(RandRef_out,mode); //w
      bgzf_mt(bgzf_fp[0],mt_cores,bgzf_buf); //
      
      //generating kstring for random reference
      char INDEL_INFO[512];
      kstring_t *RandRef;
      RandRef =(kstring_t*) calloc(1,sizeof(kstring_t));
      
      for (int i = 0; i < ChrNo; i++){
        RandRef->s = NULL;
        RandRef->l = RandRef->m = 0;
        ksprintf(RandRef,">RandChr%d\n",i+1);
        rand_val = mrand_pop_long(mr);
        for (size_t l = 1; l < RandRefLen; l++){
          ksprintf(RandRef,"%c",bases[(int)(mrand_pop_long(mr) % ntval)]);
          //fprintf(stderr,"Length %lu \t seed %d \t chrNo works %d \t nt %d \t pos file %s \t Ref file %s \t base %c\n",RandRefLen,seed,i,DoNnt,RandRef_out,bases[(int)(mrand_pop_long(mr) %4)]);
        }
        ksprintf(RandRef,"\n");
        assert(bgzf_write(bgzf_fp[0],RandRef->s,RandRef->l)!=0);
      }

      free(RandRef->s);
      free(RandRef);
      bgzf_close(bgzf_fp[0]);
      free(bgzf_fp);
    }  

    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    Reference_Variant(argc,argv); 
    return 0;
}
#endif

/*
g++ RandRef.cpp ../mrand.o ../fasta_sampler.o ../RandSampling.o -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto /projects/lundbeck/scratch/wql443/htslib/libhts.a -D __WITH_MAIN__ -o RandRef

./RandRef -l 1000000 -s 1 --chrNo 1 -nn -o /projects/lundbeck/scratch/wql443/BiasProject/Simulation/RandRefChr1S1.fa
*/