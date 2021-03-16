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

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>          // std::mutex mtx;

std::mutex mtx; 

void chrdata(faidx_t *seq_ref ,int chr_no,FILE *fp){
  const char *name = faidx_iseq(seq_ref,chr_no);
  int name_len =  faidx_seq_len(seq_ref,name);
    
  //fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
  mtx.lock();
  char *data = fai_fetch(seq_ref,name,&name_len);

  int start_pos = 1;
  int end_pos = name_len; //30001000
  char seqmod[1024] = {0}; //setting it to zero it for some reason doesnt fuck up with '?' marks at the end of each try

  while(start_pos <= end_pos){
    int readlength = 10;
    int stop = start_pos + readlength;
    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);
      
    fprintf(fp,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
    fprintf(fp,"%s\n",seqmod);
    
    start_pos += readlength + 1;
  }
  mtx.unlock();
}


void files(const char* fastafile) 
{
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;

  std::cout << faidx_nseq(seq_ref) << std::endl;
  for (int i = 0; i < faidx_nseq(seq_ref); i++){
    char buffer [50];
    int n;
    sprintf(buffer,"threads/filetest%d.fa", i);
    FILE *fp;
    fp = fopen(buffer,"wb");
    std::thread first(chrdata,seq_ref,i,fp);
    first.join();
    fclose(fp);
  }
}

void foo(const char* fastafile,FILE *fp) 
{
  std::vector<std::thread> threads;

  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;

  std::thread first(chrdata,seq_ref,0,fp);
  first.join();
  std::thread second(chrdata,seq_ref,1,fp);
  second.join();
  std::thread third(chrdata,seq_ref,2,fp);
  third.join();
  std::thread fourth(chrdata,seq_ref,3,fp);
  fourth.join();
  std::thread fifth(chrdata,seq_ref,4,fp);
  fifth.join();
  std::thread sixth(chrdata,seq_ref,5,fp);
  sixth.join();
  std::thread seventh(chrdata,seq_ref,6,fp);
  seventh.join();
  std::thread eigth(chrdata,seq_ref,7,fp);
  eigth.join();
  std::thread ninth(chrdata,seq_ref,8,fp);
  ninth.join();
  //7.61 s
}

void foo2(const char* fastafile,FILE *fp) 
{
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;

  std::cout << faidx_nseq(seq_ref) << std::endl;
  for (int i = 0; i < faidx_nseq(seq_ref); i++){
    std::thread first(chrdata,seq_ref,i,fp);
    first.join();
  }
}

void foo3(const char* fastafile,FILE *fp,int thread_no) 
{
  std::vector<std::thread> threads;

  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;
  for (int i = 0; i < 15; i++){threads.push_back(std::thread(chrdata,seq_ref,i,fp));}
  for (auto &th : threads) { th.join(); }
  threads.clear();
}

void foo4(const char* fastafile,FILE *fp,int thread_no) 
{
  std::vector<std::thread> threads;

  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  
  int chr_total = 9;//faidx_nseq(seq_ref);
  int chr_idx = 0;
  int chr_no = 0;

  while(chr_idx < chr_total){
    for (chr_no; chr_no < std::min(chr_idx+thread_no,chr_total);chr_no++){
      threads.push_back(std::thread(chrdata,seq_ref,chr_no,fp));
      std::cout << "chromosome elements "<< chr_no << std::endl;
    }
    std::cout << "end of for loop" << std::endl;
    for (auto &th : threads) { th.join(); }
    threads.clear();
    chr_idx = chr_idx + thread_no;
  }
  std::cout << "end of while lopp" << std::endl;
  //for (auto &th : threads) { th.join(); }
  //threads.clear();
}

void foo5(const char* fastafile,FILE *fp,int thread_no,std::vector<std::thread> *threads) 
{

  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  
  int chr_total = 9;//faidx_nseq(seq_ref);
  int chr_idx = 0;
  int chr_no = 0;

  while(chr_idx < chr_total){
    for (chr_no; chr_no < std::min(chr_idx+thread_no,chr_total);chr_no++){
      threads -> push_back(std::thread(chrdata,seq_ref,chr_no,fp));
      std::cout << "chromosome elements "<< chr_no << std::endl;
    }
    std::cout << "end of for loop" << std::endl;
    chr_idx = chr_idx + thread_no;
  }
  std::cout << "end of while lopp" << std::endl;
}

//cat test.fa | grep '>' | cut -f1 -d: | sort -u

int main(){
  clock_t tStart = clock();
  FILE *fp;
  fp = fopen("threads/test.fa","wb");
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  //const char *fastafile = "test.fa";
  //foo(fastafile,fp);
  //foo2(fastafile,fp);
  //foo(fastafile,fp);
  //files(fastafile);
  foo4(fastafile,fp,5);
  
  fclose(fp);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}

// g++ Threads.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// cat test.fa | grep '>' | cut -c 1-6 | sort -u