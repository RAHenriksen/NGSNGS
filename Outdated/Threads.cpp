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

pthread_mutex_t data_mutex;

struct pars_for_a_thread{
  int chr_idx;
  kstring_t *fqresult;
  faidx_t *seq_ref;
};

void* thread_runner(void *arg){
  pars_for_a_thread *struct_obj = (pars_for_a_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);
  
  int start_pos = 1;
  int end_pos = chr_len; //30001000
  char seqmod[1024] = {0};

  while(start_pos <= end_pos){
    int readlength = 100;
    int stop = start_pos + readlength;
    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);
    ksprintf(struct_obj->fqresult,"@%s:%d-%d_length:%d\n%s\n+\nIIIII\n",chr_name,start_pos,start_pos+readlength,readlength,seqmod);    
    start_pos += readlength + 1;
  }

  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  
  return NULL;
}

int main(){

  //Loading in an creating my objects for the sequence files.
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  
  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);

  //creating an array with the arguments to create multiple threads;
  int nthreads=chr_no;
  pthread_t mythreads[nthreads];
  pars_for_a_thread struct_for_threads[nthreads];

  //initialzie values that should be used for each thread
  for(int i=0;i<nthreads;i++){
    struct_for_threads[i].chr_idx = i;
    struct_for_threads[i].seq_ref = seq_ref;
    struct_for_threads[i].fqresult=new kstring_t;
    struct_for_threads[i].fqresult -> l = 0;
    struct_for_threads[i].fqresult -> m = 0;
    struct_for_threads[i].fqresult -> s = NULL;
  }

  //launch all worker threads
  for(int i=0;i<nthreads;i++){
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    pthread_create(&mythreads[i],&attr,thread_runner,&struct_for_threads[i]);
  }
    
    //right now 22 thread are running at the same time
    //now join, this means waiting for all worker threads to finish
  for(int i=0;i<nthreads;i++){
    pthread_join(mythreads[i],NULL);
  }

  pthread_mutex_destroy(&data_mutex);    
    //now all are done and we only write in main program.
    //now write everything
  FILE *fp;
  fp = fopen("threads_test.fa","wb");
  for(int i=0;i<nthreads;i++){
    fprintf(fp,"%s",struct_for_threads[i].fqresult->s);
  }
  fclose(fp);
}

// g++ TK_threads.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// cat test.fa | grep '>' | cut -c 1-6 | sort -u
