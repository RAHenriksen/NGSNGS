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
  size_t ModulusNo;
  const char* PosFile;
  const char *RefVar_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Stochastic variation on input reference genome\n");
  fprintf(fp,"Usage\n./RandVar -i <Reference> -s <seed> -n <Number of random variations> -p <Position information> -o <Reference genome with variations>\n");
  fprintf(fp,"\nExample\n./RandVar -i RandRefChr1S1.fa -s 10 -n 1000000 -p position.txt -o RandRefChr1S1Var.fa\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h | --help: \t Print help page.\n");
  fprintf(fp,"-v | --version:  Print version.\n\n");
  fprintf(fp,"-i | --input: \t Reference genome in FASTA format.\n");
  fprintf(fp,"-s | --seed: \t Seed for random generators.\n");
  fprintf(fp,"-n | --number: \t Number of stochastic variations to be included in the input reference, conflicts with -m|--modulus.\n");
  fprintf(fp,"-m | --modulus:  Every N'th position to be altered, conflicts with -n|--number.\n");
  fprintf(fp,"-p | --pos: \t Output chromomsomal coordinates for the incorporated variations in plain text.\n");
  fprintf(fp,"-o | --output: \t Altered reference genomes saved as .fa\n");
  exit(1);
  return 0;
}
 

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->Ref_input = NULL;
  mypars->seed = 0;
  mypars->VarNumber = 0;
  mypars->ModulusNo = 0;
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
    else if(strcasecmp("-m",*argv)==0 || strcasecmp("--modulus",*argv)==0){
      mypars->ModulusNo = atol(*(++argv));
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

// Comparison function for sorting
int SortChr(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}

// Comparison function for sorting
int SortChrPos(const void *a, const void *b) {
    const char *str_a = *(const char **)a;
    const char *str_b = *(const char **)b;

    // Extract column 1 (sequence name)
    char seq_name_a[50];
    char seq_name_b[50];
    sscanf(str_a, "%s", seq_name_a);
    sscanf(str_b, "%s", seq_name_b);

    // Compare based on column 1 (sequence name)
    int result = strcmp(seq_name_a, seq_name_b);
    if (result != 0) {
        return result;
    }

    // Extract column 2 (position)
    int value_a, value_b;
    sscanf(str_a, "%*s %d", &value_a);
    sscanf(str_b, "%*s %d", &value_b);

    // Compare based on column 2 (position)
    return value_a - value_b;
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
      size_t ModulusNo = mypars->ModulusNo;
      const char* Coord_file = mypars->PosFile;
      const char* SubsetChr = NULL;

      const char *bases = "ACGTN";
      fasta_sampler *fs = fasta_sampler_alloc(Ref_input,SubsetChr);
      //fprintf(stderr,"\t-> Chromosome name %s \t length %d\n",fs->seqs_names[fs->nref-1],fs->seqs_l[fs->nref-1]);
      mrand_t *mr = mrand_alloc(0,seed);
      
      //int chr_idx = fs->nref-1;
      int chr_idx;
      //Choose random chromosome index each time
     
      long rand_val;
      char buf[1024];

      FILE *fp;
      fp = fopen(Coord_file, "w");

      char **data = (char **)malloc(var_operations * sizeof(char *));
      int data_count = 0;

      if(var_operations > 0){
        for (int i = 0; i < var_operations;){
          chr_idx = (int)(mrand_pop_long(mr) % (fs->nref));
          rand_val = mrand_pop_long(mr);
          int pos = (int)(abs(rand_val) % fs->seqs_l[chr_idx]);
          //fprintf(stderr,"Random value %d with seed %d \t %ld \t and position %d\n",i,seed,rand_val,pos);
          if (fs->seqs[chr_idx][pos] != 'N'){  
            char previous; 
            previous = fs->seqs[chr_idx][pos];
            char altered;
            altered = bases[(int)(mrand_pop_long(mr) %4)];
            
            while(previous == altered){
              altered = bases[(int)(mrand_pop_long(mr) %4)];
            }
            fs->seqs[chr_idx][pos] = altered;

            char entry[100]; // Adjust the buffer size as needed
            // Append data to the array
            sprintf(entry, "%s \t %lu \t ref \t %c \t alt \t %c",fs->seqs_names[chr_idx], pos+1, previous, altered);

            if (data_count < var_operations) {
              data[data_count++] = strdup(entry);
            }
            //fprintf(fp,"%s \t %d \t ref \t %c \t alt \t %c\n",fs->seqs_names[chr_idx],pos,previous,altered);
            i++;
          }
          else
          {
            continue;
          }
        }
      }
      else if(ModulusNo > 0){
        for (size_t moduluspos = 0; moduluspos <= fs->seqs_l[fs->nref-1];moduluspos++){
          chr_idx = (int)(mrand_pop_long(mr) % (fs->nref));
          if (moduluspos % ModulusNo == 0){
            //fprintf(stderr,"check %lu \t %c\n",moduluspos,fs->seqs[chr_idx][moduluspos]);
            if (fs->seqs[chr_idx][moduluspos] != 'N'){
              char previous; 
              previous = fs->seqs[chr_idx][moduluspos];
              char altered;
              altered = bases[(int)(mrand_pop_long(mr) %4)];
              
              while(previous == altered){
                //fprintf(stderr,"Chromosome \t %s \t length %d \t position \t %d \t original \t %c \t altered \t %c\n",fs->seqs_names[chr_idx],fs->seqs_l[chr_idx],moduluspos,previous,altered);
                altered = bases[(int)(mrand_pop_long(mr) %4)];
                //fprintf(stderr,"Chromosome \t %s \t length %d \t position \t %d \t original \t %c \t altered \t %c\n",fs->seqs_names[chr_idx],fs->seqs_l[chr_idx],moduluspos,previous,altered);
              }
              fs->seqs[chr_idx][moduluspos] = altered;

              char entry[100]; // Adjust the buffer size as needed
              // Append data to the array
              sprintf(entry, "%s \t %lu \t ref \t %c \t alt \t %c",fs->seqs_names[chr_idx], moduluspos+1, previous, altered);

              if (data_count < var_operations) {
                data[data_count++] = strdup(entry);
              }
            }
          }
          else
          {
            continue;
          }
        }
      }
      else if(ModulusNo > 0 && var_operations > 0){
        fprintf(stderr,"Only use either -n or -m");exit(0);
      }

      // Sort the data
      qsort(data, data_count, sizeof(char *), SortChrPos);

      // Write the sorted data to the file
      for (int i = 0; i < data_count; i++) {
          fprintf(fp, "%s\n", data[i]);
          free(data[i]); // Free the dynamically allocated memory
      }
      fclose(fp);
      
      // Free the array
      free(data);

      for(int i=0;i<fs->nref;i++){
        snprintf(buf,1024,"%s_ngsngs",fs->seqs_names[i]);
        fs->seqs_names[i] = strdup(buf);
      }
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
