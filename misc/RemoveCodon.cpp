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

typedef struct{
  const char *Ref_input;
  const char *RefVar_out;
  const char *Codon;
  const char *Codonout;
  int StopCodon;
  int seed;
  int iter;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Iteratively removal of specific codons\n");
  fprintf(fp,"Usage\n./RemoveCodon -i <Reference> -s <seed> -iter <No. Replacment iterations> --codon <Codon to be replaced> --codonout <The desired codon> -o <Reference genome with removed codons>\n");
  fprintf(fp,"\nExample\n./RemoveCodon -i RandRefChr1S1Var.fa -o RandRefChr1S1Stop.fa -s 10 -iter 10 --codon AGT --codonout GAT\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h    | --help: \t Print help page.\n");
  fprintf(fp,"-v    | --version: \t Print version.\n\n");
  fprintf(fp,"-i    | --input: \t Reference genome in FASTA format.\n");
  fprintf(fp,"-o    | --output: \t Altered reference genomes in FASTA format, without the specified codons.\n");
  fprintf(fp,"--codon \t\t Specify which codon should be removed, if used alone a random replacement will occur.\n");
  fprintf(fp,"--codonout \t\t Specify which --codon should be removed, and replaced by specified --codonout.\n");
  fprintf(fp,"--stop \t\t\t Replace all stop codons with randomly generated codons.\n");
  fprintf(fp,"-s    | --seed: \t seed for random generators.\n");
  fprintf(fp,"-iter | --iterations: \t Number of times the codon replacements are performed (Due to triplets, replacement of one stop codon might lead to another), default = 5.\n");
  exit(1);
  return 0;
}
 
argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->Ref_input = NULL;
  mypars->RefVar_out = NULL;
  mypars->Codon = NULL;
  mypars->StopCodon = 0;
  mypars->seed = 0;
  mypars->iter = 5;

  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Ref_input = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->RefVar_out = strdup(*(++argv));
    }
    else if(strcasecmp("--codon",*argv)==0){
      mypars->Codon = strdup(*(++argv));
    }
    else if(strcasecmp("--codonout",*argv)==0){
      mypars->Codonout = strdup(*(++argv));
    }
    else if(strcasecmp("--stop",*argv)==0){
      mypars->StopCodon = 1;
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->seed = atoi(*(++argv));
    }
    else if(strcasecmp("-iter",*argv)==0 || strcasecmp("--iterations",*argv)==0){
      mypars->iter = atoi(*(++argv));
    }    
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  /*
  else if(strcasecmp("-cf",*argv)==0 || strcasecmp("--codonfile",*argv)==0){
      mypars->CodonFile = strdup(*(++argv));
    }
  */
  return mypars;
}

void RemoveStop(fasta_sampler *fs,int seed){
  // Function to Replace stop codons with another random codon.
  
  for(int chr_no = 0; chr_no<fs->nref;chr_no++){
    mrand_t *mr = mrand_alloc(0,seed+chr_no);
    char* tagPtr = strstr(fs->seqs[chr_no], "TAG");
    char* taaPtr = strstr(fs->seqs[chr_no], "TAA");
    char* tgaPtr = strstr(fs->seqs[chr_no], "TGA");

    long rand_val;
    rand_val = mrand_pop_long(mr);

    while (tagPtr != NULL) {
      int replacementChoice = (int)(mrand_pop_long(mr) % 2); // Generate a random number between 0 and 2
      switch (replacementChoice){
        case 0:
          strncpy(tagPtr, "TAC", 3);
          break;
        case 1:
          strncpy(tagPtr, "TAT", 3);
          break;
      }
      tagPtr = strstr(fs->seqs[chr_no], "TAG");
    }

    while (taaPtr != NULL) {
      int replacementChoice = (int)(mrand_pop_long(mr) % 2); // Generate a random number between 0 and 2
      switch (replacementChoice){
        case 0:
          strncpy(taaPtr, "TAC", 3);
          break;
        case 1:
          strncpy(taaPtr, "TAT", 3);
          break;
      }
      taaPtr = strstr(fs->seqs[chr_no], "TAA");
    }

    while (tgaPtr != NULL) {
      int replacementChoice = (int)(mrand_pop_long(mr) % 3); // Generate a random number between 0 and 2
      switch (replacementChoice){
        case 0:
          strncpy(tgaPtr, "TGC", 3);
          break;
        case 1:
          strncpy(tgaPtr, "TGG", 3);
          break;
        case 2:
          strncpy(tgaPtr, "TGT", 3);
          break;
      }
      tgaPtr = strstr(fs->seqs[chr_no], "TGA");
    }
  }
}

void RemoveCodon(fasta_sampler *fs,int seed,const char* codon){
  // Function to Replace stop codons with another random codon.
  
  for(int chr_no = 0; chr_no<fs->nref;chr_no++){
    mrand_t *mr = mrand_alloc(0,seed+chr_no);
    char* codonPtr = strstr(fs->seqs[chr_no], codon);

    long rand_val;
    rand_val = mrand_pop_long(mr);

    char newString[4];  // Allocate space for the new string (including null terminator)

    newString[0] = codon[0]; 
    newString[1] = codon[1];

    const char *myString = newString; 

    while (codonPtr != NULL){
      if(codon[2] == 'A'){
        const char *bases = "CGT";
        newString[2] = bases[(int)(mrand_pop_long(mr) % 3)];
        newString[3] = '\0';
        myString = newString;
        strncpy(codonPtr,myString, 3);
        codonPtr = strstr(fs->seqs[chr_no], codon);
      }
      else if(codon[2] == 'G'){
        const char *bases = "ACT";
        newString[2] = bases[(int)(mrand_pop_long(mr) % 3)];
        newString[3] = '\0';
        myString = newString;
        strncpy(codonPtr,myString, 3);
        codonPtr = strstr(fs->seqs[chr_no], codon);
      }
      else if(codon[2] == 'C'){
        const char *bases = "AGT";
        newString[2] = bases[(int)(mrand_pop_long(mr) % 3)];
        newString[3] = '\0';
        myString = newString;
        strncpy(codonPtr,myString, 3);
        codonPtr = strstr(fs->seqs[chr_no], codon);
      }
      else if(codon[2] == 'T'){
        const char *bases = "ACG";
        newString[2] = bases[(int)(mrand_pop_long(mr) % 3)];
        newString[3] = '\0';
        myString = newString;
        strncpy(codonPtr,myString, 3);
        codonPtr = strstr(fs->seqs[chr_no], codon);
      }
    }
  }
}

void RemoveSpecificCodon(fasta_sampler *fs,int seed,const char* codon,const char* codonout){
  // Function to Replace stop codons with another random codon.
  
  for(int chr_no = 0; chr_no<fs->nref;chr_no++){
    mrand_t *mr = mrand_alloc(0,seed+chr_no);
    char* codonPtr = strstr(fs->seqs[chr_no], codon);
    while (codonPtr != NULL){
      strncpy(codonPtr,codonout, 3);
      codonPtr = strstr(fs->seqs[chr_no], codon);
    }
  }
}

int ProcessRef(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    else{
      mypars = getpars(argc,argv);
      const char* Ref_input = mypars->Ref_input;
      const char* RefVar_out = mypars->RefVar_out;
      const char* Codon = mypars->Codon;
      const char* Codonout = mypars->Codonout;
      int StopCodon = mypars->StopCodon;
      //const char *CodonFile = mypars->CodonFile;
      int seed = mypars->seed;
      int iter = mypars->iter;

      const char *bases = "ACGTN";

      const char* SubsetChr = NULL;

      fasta_sampler *fs = fasta_sampler_alloc(Ref_input,SubsetChr);

      // name fs->seqs_names[fs->nref-1]
      // length of chr fs->seqs_l[fs->nref-1]
      // sequence of chr fs->seqs_l[chr_idx]

      if(StopCodon > 0){
        fprintf(stderr,"Replacing all stop codons with random codons not belonging to the stop codons\n");
        for(int i=0;i<iter;i++){
          RemoveStop(fs,seed+i);
        }
      }
      else if(Codon != NULL && Codonout == NULL){
        fprintf(stderr,"Replacing specified codon with a random codon\n");
        for(int i=0;i<iter;i++){
          //fprintf(stderr,"seed %d \t it %d \t iter %d\n",seed,i,iter);
          RemoveCodon(fs,seed+i,Codon);
        }
      }
      else if(Codon != NULL && Codonout != NULL){
        fprintf(stderr,"Replacing specified codon with another specified codon\n");
        for(int i=0;i<iter;i++){
          RemoveSpecificCodon(fs,seed+i,Codon,Codonout);
        }
      }
      else if((StopCodon > 0 && Codon != NULL) || (StopCodon > 0 && Codonout != NULL)){
        fprintf(stderr,"Only provide either --stop or --codon");
        exit(1);
      }
      else if((Codon == NULL && Codonout != NULL)){
        fprintf(stderr,"When providing --codonout also provide --codon");
        exit(1);
      }

      char buf[1024];
      for(int i=0;i<fs->nref;i++){
        snprintf(buf,1024,"%s_rm",fs->seqs_names[i]);
        fs->seqs_names[i] = strdup(buf);
      }
      dump_internal(fs,RefVar_out);
      
      fasta_sampler_destroy(fs);
    }  
    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    ProcessRef(argc,argv); 
    return 0;
}
#endif

/*
g++ RemoveCodon.cpp ../mrand.o ../fasta_sampler.o ../RandSampling.o -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto /willerslev/users-shared/science-snm-willerslev-wql443/WP1/htslib/libhts.a -D __WITH_MAIN__ -o ChrRandVar

./ChrRandVar -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -s 10 -n 1000000 -p test.txt -o output.fa
*/
