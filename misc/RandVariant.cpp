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
  int seed;
  size_t VarNumber;
  const char* PosFile;
  const char *RefVar_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Stochastic variation on input reference genome\n");
  fprintf(fp,"Usage\n./ChrRandVar -i <Reference> -s <seed> -n <Number of random variations> -p <Position information> -o <Reference genome with variations>\n");
  fprintf(fp,"\nExample\n./ChrRandVar -i chr22.fa -s 10 -n 1000000 -p position.txt -o output.fa\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-i   | --input: \t\t Reference genome .fa format\n");
  fprintf(fp,"-s   | --seed: \t\t\t seed for random generators\n");
  fprintf(fp,"-n   | --number: \t\t Number of stochastic variations to be included in the input reference\n");
  fprintf(fp,"-p   | --pos: \t\t\t Internal txt file with chromomsomal coordinate for the incorporated variations\n");
  fprintf(fp,"-o   | --output: \t\t Altered reference genomes saved as .fa\n");
  exit(1);
  return 0;
}
 

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->Ref_input = NULL;
  mypars->seed = 0;
  mypars->VarNumber = 0;
  mypars->PosFile = NULL;
  mypars->RefVar_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Ref_input = strdup(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->seed = atoi(*(++argv));
    }
    else if(strcasecmp("-n",*argv)==0 || strcasecmp("--number",*argv)==0){
      mypars->VarNumber = atol(*(++argv));
    }
    else if(strcasecmp("-p",*argv)==0 || strcasecmp("--pos",*argv)==0){
      mypars->PosFile = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->RefVar_out = strdup(*(++argv));
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
      const char* Ref_input = mypars->Ref_input;
      const char* RefVar_out = mypars->RefVar_out;
      int seed = mypars->seed;
      size_t var_operations = mypars->VarNumber;
      const char* Coord_file = mypars->PosFile;
      const char* SubsetChr = NULL;
      
      const char *bases = "ACGTN";
      fasta_sampler *fs = fasta_sampler_alloc(Ref_input,SubsetChr);
      fprintf(stderr,"\t-> Chromosome name %s \t length %d\n",fs->seqs_names[fs->nref-1],fs->seqs_l[fs->nref-1]);
      mrand_t *mr = mrand_alloc(0,seed);
      int chr_idx = fs->nref-1;
      long rand_val;

      char buf[1024];

      FILE *fp;
	    fp = fopen(Coord_file, "w");

      for (int i = 0; i < var_operations;){
        rand_val = mrand_pop_long(mr);
        int pos = (int)(abs(rand_val) % fs->seqs_l[chr_idx]);
        //fprintf(stderr,"Random value %d with seed %d \t %ld \t and position %d\n",i,seed,rand_val,pos);
        if (fs->seqs[chr_idx][pos] != 'N'){  
          char previous = fs->seqs[chr_idx][pos];
          fs->seqs[chr_idx][pos] = bases[(int)(mrand_pop_long(mr) %4)];
          char altered = fs->seqs[chr_idx][pos];
          fprintf(fp,"Chromosome \t %s \t length %d \t position \t %d \t original \t %c \t altered \t %c\n",fs->seqs_names[chr_idx],fs->seqs_l[chr_idx],pos,previous,altered);
          i++;
        }
        else
        {
          continue;
        }
      }
      fclose(fp);
      snprintf(buf,1024,"%s_ngsngs",fs->seqs_names[chr_idx]);
      fs->seqs_names[fs->nref-1] = strdup(buf);
      dump_internal(fs,RefVar_out);

      fasta_sampler_destroy(fs);
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
g++ RandVariant.cpp ../mrand.o ../fasta_sampler.o ../RandSampling.o -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto /willerslev/users-shared/science-snm-willerslev-wql443/WP1/htslib/libhts.a -D __WITH_MAIN__ -o ChrRandVar

./ChrRandVar -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -s 10 -n 1000000 -p test.txt -o output.fa
*/