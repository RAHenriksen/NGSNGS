#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <cstdlib>
#include <ctime>

#include <cstdio>
#include <cassert>
#include <cstdint>

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>          // std::mutex mtx;

#include "SimulAncient_func.h"

pthread_mutex_t data_mutex;

// -------------- SINGLE END DATA ------------- //

struct Parsarg_for_Fafq_se_thread{
  int chr_idx;
  kstring_t *fqresult_r1;
  faidx_t *seq_ref;
  const char* Ill_err;
  const char* read_err_1;
  bool Adapter_flag;
  const char* Adapter_1;
};

void* Fafq_thread_se_run(void *arg){
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);

  //just load in first file to count number of lines
  std::ifstream file(struct_obj->read_err_1);
  int Line_no = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();
  size_t lib_size = Line_no/4; 

  // loading in error profiles for nt substitions and read qual creating 2D arrays
  double** Error_2darray = create2DArray(struct_obj->Ill_err,4,280);
  double** R1_2Darray = create2DArray(struct_obj->read_err_1,8,Line_no);

  // creating random objects for all distributions.
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);

  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  char seqmod[1024] = {0};

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  
  char qual[1024] = "";
  
  // for the coverage examples
  const int Number = 100000000;
  int frag_len[Number] = {0};
  
  for (size_t i = 0; i < Number; i++){std::cout << frag_len[i] << ' ';}
  int cov = 1;
  float cov_current = 0;
  int start_pos;
  srand(time(NULL));
  int length;
  
  int nread = 0;

  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Create_se_threads(const char *fastafile,faidx_t *seq_ref,const char *output,bool Adapt_flag,const char* Adapter_1,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  int chr_no = 0;

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

  int i = 0; int j = 0; int k = 0;
  while(chr_idx < nthreads){
    //initialzie values that should be used for each thread
    for (i; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      std::cout << "chr_no" << i << std::endl;
    }
    
    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafq_thread_se_run,&struct_for_threads[j]);
    }

    //now join, this means waiting for all worker threads to finish
    for(k;k<std::min(chr_idx+thread_no,nthreads);k++){
      pthread_join(mythreads[k],NULL);
    }
    std::cout << "end of for loop" << std::endl;
    chr_idx = chr_idx + thread_no;
  }

  pthread_mutex_destroy(&data_mutex);    
    //now all are done and we only write in main program.
    //now write everything

  if (std::strcmp(output, "gz") == 0){
    gzFile gz1;
    gz1 = gzopen("threads_test_cov.fq.gz","wb");
    for(int i=0;i<nthreads;i++){
      gzwrite(gz1,struct_for_threads[i].fqresult_r1->s,struct_for_threads[i].fqresult_r1->l);
      struct_for_threads[i].fqresult_r1->l = 0;
    }
    gzclose(gz1);
  }
  else if (std::strcmp(output, "fq") == 0){
    FILE *fp1;
    fp1 = fopen("threads_test_cov.fq","wb");
    for(int i=0;i<nthreads;i++){
      fprintf(fp1,"%s",struct_for_threads[i].fqresult_r1->s);
    }
    fclose(fp1);
  }

  return NULL;
}

// ------------------------------ //

int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa  hg19canon.fa
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  
  bool Adapt_flag;

  if(argc>4){
    int thread_to_run = atoi(argv[4]);
    if (thread_to_run > chr_total)
    {
      fprintf(stderr,"The number of threads (%i) cannot exceed the number of chromosomes (%i) in the input fasta file\n",thread_to_run,chr_total);
    }
    else
    {
      if (std::strcmp(argv[2], "SE") == 0 || std::strcmp(argv[2], "single") == 0 || std::strcmp(argv[2], "single-end") == 0)
      {
        if (std::strcmp(argv[3], "False") == 0 || std::strcmp(argv[3], "false") == 0 || std::strcmp(argv[3], "F") == 0)
        {
          Adapt_flag = false;
          const char* Adapter_1 = NULL;
          Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,thread_to_run,chr_total);
        }
        else if (std::strcmp(argv[3], "True") == 0 || std::strcmp(argv[3], "true") == 0 || std::strcmp(argv[3], "T") == 0)
        {
          Adapt_flag = true;
          const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
          Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,thread_to_run,chr_total);
        }
      }
      //Create_threads(fastafile,seq_ref,argv[1],thread_to_run,chr_total);
    }
  }
  else{
    if (std::strcmp(argv[2], "SE") == 0 || std::strcmp(argv[2], "single") == 0 || std::strcmp(argv[2], "single-end") == 0){
      if (std::strcmp(argv[3], "False") == 0 || std::strcmp(argv[3], "false") == 0 || std::strcmp(argv[3], "F") == 0){
        Adapt_flag = false;
        const char* Adapter_1 = NULL;
        Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,chr_total,chr_total);
      }
      else if (std::strcmp(argv[3], "True") == 0 || std::strcmp(argv[3], "true") == 0 || std::strcmp(argv[3], "T") == 0){
        Adapt_flag = true;
        const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
        Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,chr_total,chr_total);
      }
    }
    //Create_threads(fastafile,seq_ref,argv[1],chr_total,chr_total);
  } 
}

// g++ SimulAncient_func.cpp FaFq_thread_Cov.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
// ./a.out fq SE F 3
