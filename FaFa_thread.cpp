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
#include <bits/stdc++.h>

#include "SimulAncient_func.h"

pthread_mutex_t data_mutex;

struct Parsarg_for_Fafa_thread{
  int chr_idx;
  kstring_t *fqresult;
  faidx_t *seq_ref;
  //const char* Output_type;
};

//-------------- FUNCTIONS ----------------//

void* Fafa_thread_run(void *arg){
  Parsarg_for_Fafa_thread *struct_obj = (Parsarg_for_Fafa_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);
  
  int start_pos = 1;
  int end_pos = chr_len; //30001000
  char seqmod[1024] = {0};

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  while(start_pos <= end_pos){
    //int readlength = drand48()*(80.0-30.0)+30.0;
    int readlength = sizearray[SizeDist[1](gen)];
    int stop = start_pos+(int) readlength;

    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);
    
    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += readlength + 1;}
    else{
      char nt[] = "tT";
      Deamin_char(seqmod,nt,readlength);
      ksprintf(struct_obj->fqresult,">%s:%d-%d_length:%d\n%s\n",chr_name,start_pos+1,start_pos+readlength,readlength,seqmod);
    }

    start_pos += readlength + 1;
    memset(seqmod, 0, sizeof seqmod);
  }

  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  
  return NULL;
}

void* Create_threads(const char *fastafile,faidx_t *seq_ref,const char *output,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  int chr_no = 0;
  
  
  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  //int nthreads=atai(argv[2]);
  int nthreads = chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafa_thread struct_for_threads[nthreads];

  int i = 0; int j = 0; int k = 0;
  while(chr_idx < nthreads){
    //initialzie values that should be used for each thread
    for (i; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].fqresult=new kstring_t;
      struct_for_threads[i].fqresult -> l = 0;
      struct_for_threads[i].fqresult -> m = 0;
      struct_for_threads[i].fqresult -> s = NULL;
      std::cout << "chr_no" << i << std::endl;
    }

    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafa_thread_run,&struct_for_threads[j]);
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
    gzFile gz;
    gz = gzopen("threads_test.fa.gz","wb");
    for(int i=0;i<nthreads;i++){
      gzwrite(gz,struct_for_threads[i].fqresult->s,struct_for_threads[i].fqresult->l);
    }
    gzclose(gz);
  }
  else if (std::strcmp(output, "fa") == 0){
    FILE *fp;
    fp = fopen("threads_test.fa","wb");
    for(int i=0;i<nthreads;i++){
      fprintf(fp,"%s",struct_for_threads[i].fqresult->s);
    }
    fclose(fp);
  }

  return NULL;
}

// --------------------- //
int main(int argc,char **argv){
  // "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr1_2.fa" // hg19canon.fa
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr1_2.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);

  if(argc>2){
    int thread_to_run = atoi(argv[2]);
    if (thread_to_run > chr_total)
    {
      fprintf(stderr,"The number of threads (%i) cannot exceed the number of chromosomes (%i) in the input fasta file\n",thread_to_run,chr_total);
    }
    else
    {
      Create_threads(fastafile,seq_ref,argv[1],thread_to_run,chr_total);
    }
  }
  else
  {
    Create_threads(fastafile,seq_ref,argv[1],chr_total,chr_total);
  } 
}

// g++ SimulAncient_func.cpp FaFa_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
// ./a.out fa 5
