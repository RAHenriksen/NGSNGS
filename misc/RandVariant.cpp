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
  int generation;
  double Mutation;
  size_t VarNumber;
  const char* PosFile;
  const char *RefVar_out;
}argStruct;

//functions returns head of chromosome, with posB, posE, chr_idx and fraglength set accordingly
char *sample_full_chr(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx){
  chr_idx = 0;
  if(fs->nref>1)
    chr_idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[chr_idx];

  char *ret =fs->seqs[chr_idx]; 
  chr_idx = fs->realnameidx[chr_idx];
  return ret;
}

int HelpPage(FILE *fp){
  fprintf(fp,"Stochastic variation on input reference genome\n");
  fprintf(fp,"Usage\n./RandVar --input <Reference> --seed <seed> --mutationrate/--number <Mutation rate or number of random variations> --position <Storing the information of altered position in text based format> --output <Saving altered fasta file>\n");
  fprintf(fp,"\nExample\n./RandVar --input hg19.fa --seed 10 --mutationrate 1.0E-8 --positions position.txt --output hg19var.fa\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h | --help: \t Print help page.\n");
  fprintf(fp,"-v | --version:  Print version.\n\n");
  fprintf(fp,"-i | --input: \t Reference genome in FASTA format.\n");
  fprintf(fp,"-s | --seed: \t Seed for random generators.\n");
  fprintf(fp,"-n | --number: \t Number of stochastic variations to be included in the input reference, conflicts with -m.\n");
  fprintf(fp,"-m | --mutationrate: \tCalculate the number of position to be altered, based on a fixed mutation rate pr generation.\n");
  fprintf(fp,"-g | --generations: \t Specify the number of generations for which the fized mutation rate (-m) alters the genome, default = 1.\n");
  fprintf(fp,"-p | --positions: \t Output chromomsomal coordinates for the incorporated variations in plain text.\n");
  fprintf(fp,"-o | --output: \t Altered reference genome\n");
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
  mypars->Mutation = 0.0;
  mypars->generation =1;
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
    else if(strcasecmp("-m",*argv)==0 || strcasecmp("--mutationrate",*argv)==0){
      mypars->Mutation = atof(*(++argv));
    }
    else if(strcasecmp("-g",*argv)==0 || strcasecmp("--generations",*argv)==0){
      mypars->generation = atoi(*(++argv));
    }
    else if(strcasecmp("-p",*argv)==0 || strcasecmp("--positions",*argv)==0){
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
    fprintf(stderr,"INSIDE SortChr function\n");
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
      int generation = mypars->generation;
      size_t var_operations = mypars->VarNumber;
      const char* Coord_file = mypars->PosFile;
      const char* SubsetChr = NULL;
      double Mutation = mypars->Mutation;

      if(Mutation > 0.0 && var_operations > 0){
        fprintf(stderr,"Only use either -n or -m");
        exit(1);
      }

      if(generation > 1 && Mutation == 0.0){
        fprintf(stderr,"The fixed mutation rate (-m) for each generation (-g) is not specified\n");
        exit(1);
      }

      const char *bases = "ACGTN";
      fasta_sampler *fs = fasta_sampler_alloc(Ref_input,SubsetChr);
      //fprintf(stderr,"\t-> Chromosome name %s \t length %d\n",fs->seqs_names[fs->nref-1],fs->seqs_l[fs->nref-1]);
      mrand_t *mr = mrand_alloc(0,seed);
    
      //int chr_idx = fs->nref-1;
      int chr_idx;

      // Calculate the number of variations made
      size_t num_variations = 0;
      size_t genome_len = 0;
      
      int chr_idx_tmp;
      
      for(int idx = 0; idx < (int)fs->nref;idx++){
        genome_len += fs->seqs_l[idx];
        //fprintf(stderr,"Chromosome idx %d and name %s and total genome length %zu\n",idx,fs->seqs_names[idx],genome_len);
      }
      
      if (Mutation > 0){
        num_variations = (size_t) genome_len*Mutation*generation;
      }
      else{
        num_variations = var_operations;
      }
      fprintf(stderr,"\t-> Full genome length %zu with %zu variations to be included \n",genome_len,num_variations);
      //exit(1);
      long rand_val;
      char buf[1024];

      FILE *fp;
      fp = fopen(Coord_file, "w");

      // Allocate memory for the data i need to store the actual variations
      char **data = (char **)malloc(num_variations * sizeof(char *));
      int data_count = 0;
      char *chr;
      for (int i = 0; i < num_variations;){
        chr_idx = 0; //(int)(mrand_pop_long(mr) % (fs->nref));
        //Choose random chromosome index each time
        
      //char *sample_full_chr(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx){
        char *chrseq = sample_full_chr(fs,mr,&chr,chr_idx);
        //std::cout << chrseq[2] << std::endl;

        rand_val = mrand_pop_long(mr);
        int pos = (int)(abs(rand_val) % fs->seqs_l[chr_idx]);
        //fprintf(stderr,"Random value %d with seed %d \t %ld \t and position %d\n",i,seed,rand_val,pos);
        if (chrseq[pos] != 'N'){  
          char previous; 
          previous = chrseq[pos];
          char altered;
          altered = bases[(int)(mrand_pop_long(mr) %4)];
            
          while(previous == altered){
            altered = bases[(int)(mrand_pop_long(mr) %4)];
          }
          fs->seqs[chr_idx][pos] = altered;

          char *entry = (char *)malloc(10000); // Allocate memory for the entry
          if (entry == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return 1; // Or handle the allocation failure appropriately
          }

          // Append data to the array
          sprintf(entry, "%s \t %lu \t ref \t %c \t alt \t %c",fs->seqs_names[chr_idx], pos+1, previous, altered);
          //entry 1 RandChr2 	 120591 	 ref 	 C 	 alt 	 T
          //std::cout << "entry 1 " << var_operations << std::endl;
          if (data_count < num_variations) {
            data[data_count++] = strdup(entry);
          }
          //fprintf(fp,"%s \t %d \t ref \t %c \t alt \t %c\n",fs->seqs_names[chr_idx],pos,previous,altered);
          i++;
          free(entry);
        }
        else{
          continue;
        }
      }
      
      //std::cout << "before qsort "<< std::endl;
      // Sort the data
      qsort(data, data_count, sizeof(char *), SortChrPos);

      //std::cout << "after qsort "<< std::endl;

      // Write the sorted data to the file      
      for (int i = 0; i < num_variations; i++){
        fprintf(fp, "%s\n", data[i]);
        free(data[i]); // Free the dynamically allocated memory
      }
      fclose(fp);
      //std::cout << "after fclose "<< std::endl;

      for(int i=0;i<fs->nref;i++){
        //std::cout << "i nref " << std::endl;
        snprintf(buf,1024,"%s_ngsngs",fs->seqs_names[i]);
        fs->seqs_names[i] = strdup(buf);
      }
      //std::cout << "before dump internal " << std::endl;
      dump_internal(fs,RefVar_out);
      //std::cout << "after dump internal " << std::endl;
      //fasta_sampler_destroy(fs);
      //std::cout << "after fasta destroy" << std::endl;
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
