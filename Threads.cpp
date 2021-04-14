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


struct pars_for_a_thread{
  kstring_t *fqresult;
};

void* thread_runner(void *arg){
  pars_for_a_thread *struct_obj = (pars_for_a_thread*) arg;

  for(int i=0;i<100;i++){
    ksprintf(struct_obj->fqresult,"@name\nCGTGA\n+\nIIIII\n");
  }
  //now we have generated 1mio fq reads.
  return NULL;
}

int main(){
  int nthreads=22;
  //creating an array with the arguments to create multiple threads;
  pthread_t mythreads[nthreads];
  pars_for_a_thread struct_for_threads[nthreads];

  //initialzie values that should be used for each thread
  for(int i=0;i<nthreads;i++){
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
