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

//-------------- SE ----------------//

struct Parsarg_for_Fafa_se_thread{
  int chr_idx;
  kstring_t *fqresult;
  faidx_t *seq_ref;
  const char* Ill_err;
  //const char* Output_type;
};

void* Fafa_thread_se_run(void *arg){
  Parsarg_for_Fafa_se_thread *struct_obj = (Parsarg_for_Fafa_se_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);
  
  int start_pos = 1;
  int end_pos = chr_len; //30001000
  char seqmod[1024] = {0};

  double init = 1.0;
  double cov = 1.0;
  
  // loading in error profiles for nt substitions and read qual creating 2D arrays
  double** Error_2darray = create2DArray(struct_obj->Ill_err,4,280);
  

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::random_device rd;
  std::default_random_engine gen(rd());

  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  while(start_pos <= end_pos){
    //int readlength = drand48()*(80.0-30.0)+30.0;
    int fraglength = (int) sizearray[SizeDist[1](gen)];
    // estimates brief coverage
    int dist = init/cov * fraglength; 
    //int stop = start_pos+(int) readlength;

    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+start_pos-1,150);
    }
    // case 2
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+start_pos-1,150);
    }
    // case 3
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+start_pos-1,fraglength);
    }
    
    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += dist + 1;}
    else{
      char nt[] = "tT";
      Ill_err(seqmod,Error,gen);
      //Deamin_char(seqmod,nt,fraglength);
      ksprintf(struct_obj->fqresult,">%s:%d-%d_length:%d\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod);
    }

    start_pos += dist;
    memset(seqmod, 0, sizeof seqmod);
  }

  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  
  return NULL;
}

void* Create_se_threads(const char *fastafile,faidx_t *seq_ref,const char *output,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  //int chr_no = 0;

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  //int nthreads=atai(argv[2]);
  int nthreads = chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafa_se_thread struct_for_threads[nthreads];

  while(chr_idx < nthreads){
    //initialzie values that should be used for each thread
    for (int i = 0; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].fqresult=new kstring_t;
      struct_for_threads[i].fqresult -> l = 0;
      struct_for_threads[i].fqresult -> m = 0;
      struct_for_threads[i].fqresult -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      std::cout << "chr_no" << i << std::endl;
    }

    //launch worker threads
    for (int j = 0; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafa_thread_se_run,&struct_for_threads[j]);
    }

    //now join, this means waiting for all worker threads to finish
    for(int k = 0; k<std::min(chr_idx+thread_no,nthreads);k++){
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

//-------------- PE ----------------//

struct Parsarg_for_Fafa_pe_thread{
  int chr_idx;
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  faidx_t *seq_ref;
  const char* Ill_err;
};

void* Fafa_thread_pe_run(void *arg){
  Parsarg_for_Fafa_pe_thread *struct_obj = (Parsarg_for_Fafa_pe_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);
  
  int start_pos = 1;
  int end_pos = chr_len; //30001000
  char seqmod[1024] = {0};
  char seqmod2[1024] = {0};

  double init = 1.0;
  double cov = 1.0;
  
  // loading in error profiles for nt substitions and read qual creating 2D arrays
  double** Error_2darray = create2DArray(struct_obj->Ill_err,4,280);
  

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::random_device rd;
  std::default_random_engine gen(rd());

  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);
  
  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  while(start_pos <= end_pos){
    //int readlength = drand48()*(80.0-30.0)+30.0;
    int fraglength = (int) sizearray[SizeDist[1](gen)];
    // estimates brief coverage
    int dist = init/cov * fraglength; 
    //int stop = start_pos+(int) readlength;

    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+start_pos-1,150);
      strncpy(seqmod2,data+fraglength-1-150,150);
    }
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+start_pos-1,150);
      strncpy(seqmod2,data+start_pos+fraglength-1-150,150);
    }
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+start_pos-1,fraglength);
      strncpy(seqmod2,data+start_pos-1,fraglength);
    }
    
    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += dist + 1;}
    else{
      char nt[] = "tT";
      Ill_err(seqmod,Error,gen);
      std::reverse(seqmod2, seqmod2 + strlen(seqmod2)); // jeg bliver vel nødt til med det samme at revertere den!
      Ill_err(seqmod2,Error,gen);
      DNA_complement(seqmod2);

      //Deamin_char(seqmod,nt,fraglength);
      ksprintf(struct_obj->fqresult_r1,">%s:%d-%d_length:%d\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod);
      ksprintf(struct_obj->fqresult_r2,">%s:%d-%d_length:%d\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod);
    }

    start_pos += dist;
    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
  }

  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  
  return NULL;
}

void* Create_pe_threads(const char *fastafile,faidx_t *seq_ref,const char *output,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  //int chr_no = 0;

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  //int nthreads=atai(argv[2]);
  int nthreads = chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafa_pe_thread struct_for_threads[nthreads];

  while(chr_idx < nthreads){
    //initialzie values that should be used for each thread
    for (int i = 0; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;
      struct_for_threads[i].fqresult_r2 =new kstring_t;
      struct_for_threads[i].fqresult_r2 -> l = 0;
      struct_for_threads[i].fqresult_r2 -> m = 0;
      struct_for_threads[i].fqresult_r2 -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      std::cout << "chr_no" << i << std::endl;
    }

    //launch worker threads
    for (int j = 0; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafa_thread_pe_run,&struct_for_threads[j]);
    }

    //now join, this means waiting for all worker threads to finish
    for(int k = 0; k<std::min(chr_idx+thread_no,nthreads);k++){
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
    gzFile gz2;
    gz1 = gzopen("threads_test_R1.fa.gz","wb");
    gz2 = gzopen("threads_test_R2.fa.gz","wb");
    
    for(int i=0;i<nthreads;i++){
      gzwrite(gz1,struct_for_threads[i].fqresult_r1->s,struct_for_threads[i].fqresult_r1->l);
      gzwrite(gz2,struct_for_threads[i].fqresult_r2->s,struct_for_threads[i].fqresult_r2->l);
    }
    gzclose(gz1);
    gzclose(gz2);
  }
  else if (std::strcmp(output, "fa") == 0){
    FILE *fp1;
    FILE *fp2;
    fp1 = fopen("threads_test_R1.fa","wb");
    fp2 = fopen("threads_test_R2.fa","wb");

    for(int i=0;i<nthreads;i++){
      fprintf(fp1,"%s",struct_for_threads[i].fqresult_r1->s);
      fprintf(fp2,"%s",struct_for_threads[i].fqresult_r2->s);
    }
    fclose(fp1);
    fclose(fp2);
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
      //Create_se_threads(fastafile,seq_ref,argv[1],thread_to_run,chr_total);
      Create_pe_threads(fastafile,seq_ref,argv[1],thread_to_run,chr_total);
    }
  }
  else
  {
    //Create_se_threads(fastafile,seq_ref,argv[1],chr_total,chr_total);
    Create_pe_threads(fastafile,seq_ref,argv[1],chr_total,chr_total);
  } 
}

// g++ SimulAncient_func.cpp FaFa_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
// ./a.out fa 5
