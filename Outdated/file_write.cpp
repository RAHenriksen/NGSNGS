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
#include <bits/stdc++.h>

std::mutex mtx; 

void chrdata(faidx_t *seq_ref ,int chr_no,FILE *fp){
  const char *name = faidx_iseq(seq_ref,chr_no);
  int name_len =  faidx_seq_len(seq_ref,name);
    
  //fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
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
}

void files(const char* fastafile) {
  //This creates one thread at a time and then runs each before the next
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  
  std::cout << faidx_nseq(seq_ref) << std::endl;
  for (int i = 0; i < chr_no-1; i++){
    char buffer[50];
    int n;
    sprintf(buffer,"filetmp%d.fa", i);
    FILE *fp;
    fp = fopen(buffer,"wb");
    std::thread first(chrdata,seq_ref,i,fp);
    first.join();
    fclose(fp);
  }
}

void chrdata3(faidx_t *seq_ref ,int chr_no,FILE *fp){
  const char *name = faidx_iseq(seq_ref,chr_no);
  int name_len =  faidx_seq_len(seq_ref,name);
  
  fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
  mtx.lock();
  char *data = fai_fetch(seq_ref,name,&name_len);
  mtx.unlock();
  int start_pos = 1;
  int end_pos = name_len-1; //30001000
  char seqmod[1024] = {0}; //setting it to zero it for some reason doesnt fuck up with '?' marks at the end of each try
  //std::cout <<strncpy(seqmod,data+start_pos,20) << std::endl;
  while(start_pos <= end_pos){
    int readlength = 10;
    int stop = start_pos + readlength;
    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);
      
    fprintf(fp,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
    fprintf(fp,"%s\n",seqmod);
    
    start_pos += readlength + 1;
  }
  fclose(fp);
}

void files3(const char* fastafile) 
{
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  std::vector<std::thread> threads;
  std::cout << faidx_nseq(seq_ref) << std::endl;
  FILE *fp;
  for (int i = 0; i < chr_no-1; i++){
    char buffer [50];
    sprintf(buffer,"tmp_idx%d.fa", i);
    fp = fopen(buffer,"wb");
    threads.push_back(std::thread(chrdata3,seq_ref,i,fp));
  }
  for (auto &th : threads) { th.join(); }
}

int main(){
  clock_t tStart = clock();
  //fp = fopen("threads/test.fa","wb");
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";

  /*
  files(fastafile);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  */

  files3(fastafile);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);  
  std::string cat_str = "cat tmp_idx*.fa > Hg19_aDNA.fa";
  const char *catcmd = cat_str.c_str();
  system(catcmd);

  std::string rm_str = "rm tmp_idx*.fa";
  const char *rmcmd = rm_str.c_str();
  system(rmcmd);

  printf("Time taken2: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}

// g++ file_write.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// cat test.fa | grep '>' | cut -c 1-6 | sort -u