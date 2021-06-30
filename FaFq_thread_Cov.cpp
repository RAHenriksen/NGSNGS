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

// -------------- SINGLE END DATA ------------- //

struct Parsarg_for_Fafq_se_thread{
  int chr_idx;
  kstring_t *fqresult_r1;
  faidx_t *seq_ref;
  const char* Ill_err;
  const char* read_err_1;
};


void* Fafq_thread_se_run_1(void *arg){
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

  int start_pos = 1;
  int end_pos = chr_len; //30001000

  //coverage method2
  int* chrlencov = (int *) calloc(chr_len,sizeof(int));

  char seqmod[1024] = {0};

  // for the coverage examples
  float cov = 1.00;
  float cov_current = 0;
  int rand_start;
  int fraglength;
  int nread = 0;

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  
  char qual[1024] = "";
  //int D_i = 0;
  while (cov_current < cov) {
    int D_i = 0;
    int N = chr_len;

    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
    srand48(fraglength+std::time(nullptr));
    rand_start = lrand48() % (chr_len-fraglength-1);

    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+rand_start-1,150);
    }
    // case 2
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+rand_start-1,150);
    }
    // case 3
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+rand_start-1,fraglength);
    }

    //removes reads with NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){
      //fprintf(stderr,"pch \n");
      /*
      for (int j = rand_start; j < rand_start+fraglength; j++){chrlencov[j] += 1;}
      
      for(int i = 0; i<N ; i++){
        if (chrlencov[i] == 0){N = N-1;} //trækker -1 fra da 0 jo ikke er data -> så alle positioner med værdi forskellig for 0
        else{D_i+=chrlencov[i];} // opdaterer summen D_i = D_i+j for length
      }*/
      continue;
    } // readlength + 1
    else{
      for (int j = rand_start; j < rand_start+fraglength; j++){chrlencov[j] += 1;}
      
      // for(int i = rand_start; i<rand_start+fraglength; i++){
      for(int i = 0; i<N ; i++){
        if (chrlencov[i] == 0){N = N-1;} //trækker -1 fra da 0 jo ikke er data -> så alle positioner med værdi forskellig for 0
        else{D_i+=chrlencov[i];} // opdaterer summen D_i = D_i+j for length
      }
      Ill_err(seqmod,Error,gen);
      //std::cout << "non flag "<< std::endl;
      Read_Qual2(seqmod,qual,Qualdistr1,gen);
      ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,rand_start,rand_start+fraglength-1,fraglength,seqmod,qual);
      memset(qual, 0, sizeof(qual));  
      nread++;
      fprintf(stderr,"Number of reads %d , sum %d, count %d, coverage %f\n",nread,D_i,N,cov_current);
    }
    cov_current = (float) D_i / N;
    memset(seqmod, 0, sizeof seqmod);
    
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
  }
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Fafq_thread_se_run_2(void *arg){
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

  int start_pos = 1;
  int end_pos = chr_len; //30001000

  //coverage method2
  int* chrlencov = (int *) calloc(chr_len,sizeof(int));

  char seqmod[1024] = {0};

  // for the coverage examples
  float cov = 1.00;
  float cov_current = 0;
  int rand_start;
  int fraglength;
  int nread = 0;

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  
  char qual[1024] = "";
  //int D_i = 0;
  int size_data = 0; // between 0 and chr_len
  int D_total = 0;
  
  while (cov_current < cov) {
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
    srand48(fraglength+std::time(nullptr));
    rand_start = lrand48() % (chr_len-fraglength-1);

    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+rand_start-1,150);
    }
    // case 2
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+rand_start-1,150);
    }
    // case 3
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+rand_start-1,fraglength);
    }

    //removes reads with NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){continue;}
    else{
      for (int j = rand_start; j < rand_start+fraglength; j++){
        chrlencov[j] += 1; //adds 1 for regions, which might have reads already
        D_total += 1; // to find the total depth
        if (chrlencov[j] == 1){size_data++;} // if the value is different from 1 (more reads) then we still count that as one chromosome site with data
      }

      Ill_err(seqmod,Error,gen);
      //std::cout << "non flag "<< std::endl;
      Read_Qual2(seqmod,qual,Qualdistr1,gen);
      ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,rand_start,rand_start+fraglength-1,fraglength,seqmod,qual);
      memset(qual, 0, sizeof(qual));  
      nread++;
      fprintf(stderr,"Number of reads %d , d_total %d, size_data %d, cov_current %f\n",nread,D_total,size_data,cov_current);
    }
    cov_current = (float) D_total / (chr_len-size_data);
    memset(seqmod, 0, sizeof seqmod);
    
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
  }
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Fafq_thread_se_run_3(void *arg){
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

  int start_pos = 1;
  int end_pos = chr_len; //30001000

  char seqmod[1024] = {0};

  // for the coverage examples
  double init = 1.0;
  double cov = 1.0;

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  
  char qual[1024] = "";
  
  while(start_pos <= end_pos){
    //srand48(start_pos+std::time(nullptr)); //you could make srand48(start_pos + seed_input)
   
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    // estimates brief coverage
    int dist = init/cov * fraglength; 
        
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

    //removes reads with NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += dist+1;} // readlength + 1
    else{
      Ill_err(seqmod,Error,gen);
      //std::cout << "non flag "<< std::endl;
      Read_Qual2(seqmod,qual,Qualdistr1,gen);
      ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod,qual);
      memset(qual, 0, sizeof(qual));      
    }
    start_pos += dist;
    memset(seqmod, 0, sizeof seqmod);
  }
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Create_se_threads(const char *fastafile,faidx_t *seq_ref,int thread_no,int chr_total,int covmethod){
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
      std::cout << "chr_no" << i << std::endl;
    }
    
    //launch worker threads

    if (covmethod == 1){
      for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_create(&mythreads[j],&attr,Fafq_thread_se_run_3,&struct_for_threads[j]);
      }
    }
    else if (covmethod == 2)
    {
      for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_create(&mythreads[j],&attr,Fafq_thread_se_run_2,&struct_for_threads[j]);
      }
    }
    else if (covmethod == 3)
    {
      for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_create(&mythreads[j],&attr,Fafq_thread_se_run_3,&struct_for_threads[j]);
      }
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

  FILE *fp1;
  fp1 = fopen("Fq_cov_2.fq","wb");
  for(int i=0;i<nthreads;i++){
    fprintf(fp1,"%s",struct_for_threads[i].fqresult_r1->s);
  }
  fclose(fp1);
  return NULL;
}

// ------------------------------ //


void* full_genome_test(faidx_t *seq_ref,int chr_total){
  int genome_size = 0;
  int chr_sizes[chr_total];
  const char *chr_names[chr_total];

  for (size_t i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
  }

  std::cout << genome_size << std::endl;
  std::cout << "----------" << std::endl;
  char* genome = (char*) malloc(genome_size);
  
  for (size_t i = 0; i < chr_total; i++){

    pthread_mutex_lock(&data_mutex);
    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    pthread_mutex_unlock(&data_mutex);
    std::cout << chr_sizes[i] << std::endl;
    std::cout << chr_names[i] << std::endl;
    std::cout << strlen(data) << std::endl;
    strcat(genome,data);
  }

  char seqmod[1024] = {0};
  
  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  // creating random objects for all distributions.
  std::random_device rd;
  std::default_random_engine gen(rd()); 
  
  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  
  int rand_start;
  for (size_t i = 0; i < 10; i++)
  {
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    srand48(fraglength+std::time(nullptr));
    rand_start = lrand48() % (genome_size-fraglength-1);
    std::cout << rand_start << std::endl;
    strncpy(seqmod,genome+rand_start-1,fraglength);
    std::cout << seqmod << std::endl;
    memset(seqmod, 0, sizeof seqmod);
  }
}

// ------------------------------ //
int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa  hg19canon.fa
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  //const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  
  full_genome_test(seq_ref,chr_total);
 

  //Create_se_threads(fastafile,seq_ref,chr_total,chr_total,2);
}

//g++ SimulAncient_func.cpp FaFq_thread_Cov.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl