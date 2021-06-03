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

struct Parsarg_for_Fafq_thread{
  int chr_idx;
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  faidx_t *seq_ref;
  const char* Ill_err;
  const char* read_err_1;
  const char* read_err_2;
  bool Adapter_flag;
  const char* Platform;
  const char* Adapter_1;
  const char* Adapter_2;
};

void* Platform(int i){
  std::cout << " function " << i <<std::endl;
}

void* Fafq_thread_run(void *arg){
  Parsarg_for_Fafq_thread *struct_obj = (Parsarg_for_Fafq_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);

  if (struct_obj->Platform==(const char*)"SE"){
    Platform((int) 2);
    fprintf(stderr,"SE LORT");
  }
  std::ifstream file(struct_obj->read_err_1);
  int Line_no = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();
  size_t lib_size = Line_no/4; 

  // loading in error profiles for nt substitions and read qual creating 2D arrays
  double** Error_2darray = create2DArray(struct_obj->Ill_err,4,280);
  double** R1_2Darray = create2DArray(struct_obj->read_err_1,8,Line_no);
  double** R2_2Darray = create2DArray(struct_obj->read_err_2,8,Line_no);

  // creating random objects for all distributions.
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);
  std::discrete_distribution<> Qualdistr2[Line_no];
  Qual_dist(R2_2Darray,Qualdistr2,Line_no);
  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  int start_pos = 1;
  int end_pos = chr_len; //30001000

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
  
  while(start_pos <= end_pos){
    //std::cout << "while loop" << std::endl;
    srand(time(NULL)+start_pos);
    int readlength = drand48()*(70.0-30.0)+30.0; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    // int readlength = sizearray[SizeDist[1](gen)];
    int stop = start_pos+(int) readlength;
    
    //extracts the sequence
    strncpy(seqmod,data+start_pos-1,readlength);
    //std::cout << "sequence\t" << seqmod << std::endl;
    
    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += readlength;} // readlength + 1
    else{
      //std::cout << "else" << std::endl;
      Ill_err(seqmod,Error,gen);
      //std::cout << "error" << std::endl;
      if(struct_obj->Adapter_flag==true){        
        //std::cout << "flag" << std::endl;
        char read1N[lib_size+1];
        char read2N[lib_size+1];

        memset(read1N,'N',lib_size);
        read1N[lib_size]='\0';
        memset(read2N,'N',lib_size);
        read2N[lib_size]='\0';

        char read1[lib_size + 1];
        char read2[lib_size + 1];
        //Copies sequence into both reads
        strcpy(read1, seqmod);
        strcpy(read2, seqmod);

        //creates read 1
        std::strcat(read1,struct_obj->Adapter_1); //add adapter to read
        std::strncpy(read1N, read1, strlen(read1)); //copy read+adapter into N...

        DNA_complement(read2);
        std::reverse(read2, read2 + strlen(read2));
        std::strcat(read2,struct_obj->Adapter_2);
        std::strncpy(read2N, read2, strlen(read2)); //copy read+adapter into N...

        //Creating read quality
        Read_Qual2(read1N,qual,Qualdistr1,gen);
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength-1,readlength,read1N,qual);
        memset(qual, 0, sizeof(qual));    

        Read_Qual2(read2N,qual,Qualdistr2,gen);
        ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength-1,readlength,read2N,qual);
        memset(qual, 0, sizeof(qual));    
      }
      else{
        //std::cout << "non flag "<< std::endl;
        Read_Qual2(seqmod,qual,Qualdistr1,gen);
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength-1,readlength,seqmod,qual);
        memset(qual, 0, sizeof(qual));
        DNA_complement(seqmod);
        std::reverse(seqmod, seqmod + strlen(seqmod));
        Read_Qual2(seqmod,qual,Qualdistr2,gen);
        ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength-1,readlength,seqmod,qual);
        memset(qual, 0, sizeof(qual));        
      }
    }
    start_pos += readlength;
    memset(seqmod, 0, sizeof seqmod);
  }
  //std::cout << "out of while" << std::endl;
  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Create_pe_threads(const char *fastafile,faidx_t *seq_ref,const char *output,const char* PlatformType,bool Adapt_flag,const char* Adapter_1,const char* Adapter_2,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  int chr_no = 0;

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafq_thread struct_for_threads[nthreads];

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
      struct_for_threads[i].fqresult_r2=new kstring_t;
      struct_for_threads[i].fqresult_r2 -> l = 0;
      struct_for_threads[i].fqresult_r2 -> m = 0;
      struct_for_threads[i].fqresult_r2 -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      struct_for_threads[i].read_err_2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";
      struct_for_threads[i].Platform = PlatformType;
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].Adapter_2 = Adapter_2;
      std::cout << "chr_no" << i << std::endl;
    }
    
    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafq_thread_run,&struct_for_threads[j]);
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
    gzFile gz2;
    gz1 = gzopen("threads_test_R1.fq.gz","wb");
    gz2 = gzopen("threads_test_R2.fq.gz","wb");
    for(int i=0;i<nthreads;i++){
      gzwrite(gz1,struct_for_threads[i].fqresult_r1->s,struct_for_threads[i].fqresult_r1->l);
      struct_for_threads[i].fqresult_r1->l = 0;
      gzwrite(gz2,struct_for_threads[i].fqresult_r2->s,struct_for_threads[i].fqresult_r2->l);
      struct_for_threads[i].fqresult_r2->l = 0;
    }
    gzclose(gz1);
    gzclose(gz2);
  }
  else if (std::strcmp(output, "fq") == 0){
    FILE *fp1;
    FILE *fp2;
    fp1 = fopen("threads_test_R1.fq","wb");
    fp2 = fopen("threads_test_R2.fq","wb");
    for(int i=0;i<nthreads;i++){
      fprintf(fp1,"%s",struct_for_threads[i].fqresult_r1->s);
      fprintf(fp2,"%s",struct_for_threads[i].fqresult_r2->s);
    }
    fclose(fp1);
    fclose(fp2);
  }

  return NULL;
}



int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa  hg19canon.fa
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);

  bool Adapt_flag;
  const char* PlatformType;
  PlatformType = "SE";

  int thread_to_run = 3;
  if (std::strcmp(argv[2], "False") == 0 || std::strcmp(argv[2], "false") == 0 || std::strcmp(argv[2], "F") == 0){
    Adapt_flag = false;
    const char* Adapter_1 = NULL;
    const char* Adapter_2 = NULL;
    Create_pe_threads(fastafile,seq_ref,argv[1],PlatformType,Adapt_flag,Adapter_1,Adapter_2,thread_to_run,chr_total);
  }
  else if (std::strcmp(argv[2], "True") == 0 || std::strcmp(argv[2], "true") == 0 || std::strcmp(argv[2], "T") == 0){
    Adapt_flag = true;
    const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    const char* Adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
    Create_pe_threads(fastafile,seq_ref,argv[1],PlatformType,Adapt_flag,Adapter_1,Adapter_2,
                    thread_to_run,chr_total);
  }
  //fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
}

// g++ SimulAncient_func.cpp FaFq_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
// ./a.out fq F
