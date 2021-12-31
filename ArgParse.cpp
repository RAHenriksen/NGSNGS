#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"
#include "Sampling.h"

typedef struct{
  int threads;
  size_t reads;
  int Glob_seed;
  const char *OutFormat;
  const char *OutName;
  const char *Reference;
  const char *Seq;
  const char *Adapter1;
  const char *Adapter2;
  const char *QualProfile1;
  const char *QualProfile2;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 0.0.0 \n\n");
  fprintf(fp,"Usage: ./ngsngs <input reference> <numer of reads> <output file>\n\n");
  fprintf(fp,"Options: \n");
  fprintf(fp,"-h | --help: \t\t Print help page\n");
  fprintf(fp,"-i | --input: \t\t Reference file in .fasta format to sample reads from\n");
  fprintf(fp,"-r | --reads: \t\t Number of reads to simulate from\n");
  fprintf(fp,"-o | --output: \t\t Prefix of output file name, with default extension in fasta format (.fa)\n");
  fprintf(fp,"-f | --format: \t\t File format of the simulated outpur reads\n");
  fprintf(fp,"\t <.fa||.fasta>\t nucletide sequence \n \t <.fq||.fastq>\t nucletide sequence with corresponding quality score \n \t <.sam||bam>\t Sequence Alignment Map format\n");
  fprintf(fp,"-t | --threads: \t Number of threads to use for simulation\n");
  fprintf(fp,"-s | --seed: \t\t Random seed, default value being computer time\n");
  fprintf(fp,"-a1 | --adapter1: \t Adapter sequence to add for simulated reads (SE) or first read pair (PE)\n");
  fprintf(fp,"-a2 | --adapter2: \t Adapter sequence to add for second read pair (PE) \n");
  fprintf(fp,"-seq | --sequencing: \t Simulate single-end or paired-end\n");
  fprintf(fp,"\t <SE>\t single-end \n \t <PE>\t paired-end\n");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->threads = 1;
  mypars->reads = -1;
  mypars->Glob_seed = (int) time(NULL);
  mypars->OutFormat = "fa";
  mypars->OutName = "output";
  mypars->Seq = "SE";
  mypars->Reference = NULL;
  mypars->Adapter1 = NULL;
  mypars->Adapter2 = NULL;
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  while(*argv){
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->reads = atoi(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->OutName = strdup(*(++argv));
    }
    //optional
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->threads = atoi(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv));
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      mypars->Seq = strdup(*(++argv));
      if(strcasecmp("SE",mypars->Seq)!=0 && strcasecmp("PE",mypars->Seq)!=0){HelpPage(stderr);} 
    }
    else if(strcasecmp("-a1",*argv)==0 || strcasecmp("--adapter1",*argv)==0){
      mypars->Adapter1 = strdup(*(++argv));
    }
    else if(strcasecmp("-a2",*argv)==0 || strcasecmp("--adapter2",*argv)==0){
      mypars->Adapter2 = strdup(*(++argv));
    }
    else if(strcasecmp("-q1",*argv)==0 || strcasecmp("--quality1",*argv)==0){
      mypars->QualProfile1 = strdup(*(++argv));
    }
    else if(strcasecmp("-q2",*argv)==0 || strcasecmp("--quality2",*argv)==0){
      mypars->QualProfile2 = strdup(*(++argv));
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      mypars->OutFormat = strdup(*(++argv));
    }
    /*else{
      fprintf(stderr,"unrecognized input option, see NGSNGS help page\n\n");
      HelpPage(stderr);
      } */
    
    // -l || --length, -e1 +2 || --error1 +2 
    ++argv;
  }
  return mypars;
}

//else if(strcasecmp("-f",*argv)==0 || strcasecmp("--out",*argv)==0){}

// ------------------------------ //

int main(int argc,char **argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                        strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    mypars = getpars(argc,argv);
    clock_t t = clock();
    time_t t2 = time(NULL);

    const char *fastafile = mypars->Reference; //"/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);

    fprintf(stderr,"\t-> fasta load \n");
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int Glob_seed = mypars->Glob_seed; //(int) time(NULL);
    int threads = mypars->threads;
    size_t No_reads = mypars->reads; //1e1;

    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d with %d threads\n",Glob_seed,threads);
    fprintf(stderr,"\t-> Number of simulated reads: %zd\n",No_reads);

    const char* Adapt_flag;
    const char* Adapter_1;
    if (mypars->Adapter1 != NULL){
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      //Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    }
    else{Adapt_flag = "false";}

    const char* OutputFormat = mypars->OutFormat;

    int Thread_specific_Read = static_cast<int>(No_reads/threads);

    const char* filename = mypars->OutName; //"chr22_out";
    const char* Seq_Type = mypars->Seq;

    Create_se_threads(seq_ref,threads,Glob_seed,Thread_specific_Read,filename,Adapt_flag,Adapter_1,OutputFormat,Seq_Type);

    fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }
}


// g++ NGSNGS_func.cpp atomic_fq.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// valgrind --tool=memcheck --leak-check=full  --track-origins=yes ./a.out
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -f bam -o chr22
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 1000 -f fq -s 1 -o chr22_out
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq SE -f fq -o chr22
// ./a.out -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq PE -f fq -o chr22
//  make HTSSRC=/home/wql443/scratch/htslib/
// ./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq PE -f fq -o chr22

//g++ NGSNGS_func.cpp Sampling.cpp ArgParse.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl