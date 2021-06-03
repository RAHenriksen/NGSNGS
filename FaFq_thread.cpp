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

// -------------- PAIRED END DATA ------------- //

struct Parsarg_for_Fafq_pe_thread{
  int chr_idx;
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  faidx_t *seq_ref;
  const char* Ill_err;
  const char* read_err_1;
  const char* read_err_2;
  bool Adapter_flag;
  const char* Adapter_1;
  const char* Adapter_2;
};

void* Fafq_thread_pe_run(void *arg){
  Parsarg_for_Fafq_pe_thread *struct_obj = (Parsarg_for_Fafq_pe_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);

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
  char seqmod2[1024] = {0};

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
    //srand48(start_pos+std::time(nullptr)); //you could make srand48(start_pos + seed_input)
   
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    /*std::cout << "_-------------------" << std::endl;
    std::cout << "before" << std::endl;
    std::cout << " fragment length " << fraglength << std::endl;
    std::cout << "fails" << std::endl;*/
    
    // int readlength = sizearray[SizeDist[1](gen)];
    int stop = start_pos+(int) fraglength;
    
    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+start_pos-1,150);
      strncpy(seqmod2,data+fraglength-1-150,150); //it needs to be reversed!
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

    //removes reads with NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += fraglength+1;} // readlength + 1
    else{
      Ill_err(seqmod,Error,gen);
      std::reverse(seqmod2, seqmod2 + strlen(seqmod2)); // jeg bliver vel nødt til med det samme at revertere den!
      Ill_err(seqmod2,Error,gen); // SKULLE JEG BRUGE DEN SAMME SEQMOD HER?

      if(struct_obj->Adapter_flag==true){  
        if (fraglength < 150){ // for reads to be added with adapter

          char read1N[lib_size+1];
          char read2N[lib_size+1];
          char read1[lib_size + 1];
          char read2[lib_size + 1];

          memset(read1N,'N',lib_size);
          read1N[lib_size]='\0';
          memset(read2N,'N',lib_size);
          read2N[lib_size]='\0';

          //Copies sequence into both reads
          strcpy(read1, seqmod);
          std::strcat(read1,struct_obj->Adapter_1);

          /*strcpy(read1, struct_obj->Adapter_1); //adds adapter to the beginning of 5' before the sequence
          std::strcat(read1,seqmod);*/

          DNA_complement(seqmod2);
          
          if (strlen(read1) >= 150){
            std::strncpy(read1N, read1, 150);
            strcpy(read2, seqmod2);
            std::reverse(read2, read2 + strlen(read2));
            std::strcat(read2, struct_obj->Adapter_2);
            std::strncpy(read2N, read2, 150);
          }
          else{
            //std::cout << "else " << std::endl;
            //std::cout << "read1 v2 \n" << read1 << std::endl;
            std::strncpy(read1N, read1, strlen(read1)); //copy read+adapter into N...
            //std::cout << read1N << std::endl;
            strcpy(read2, seqmod2);
            std::reverse(read2, read2 + strlen(read2));
            std::strcat(read2, struct_obj->Adapter_2);
            std::strncpy(read2N, read2, strlen(read2));
          }
          
          //Creating read quality
          Read_Qual2(read1N,qual,Qualdistr1,gen);
          ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,read1N,qual);
          memset(qual, 0, sizeof(qual));    

          Read_Qual2(read2N,qual,Qualdistr2,gen);
          ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,read2N,qual);
          memset(qual, 0, sizeof(qual));    
        }
        else{
          //std::cout << "else loop" << std::endl;
          Read_Qual2(seqmod,qual,Qualdistr1,gen);
          ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod,qual);
          memset(qual, 0, sizeof(qual));
          DNA_complement(seqmod2);
          std::reverse(seqmod2, seqmod2 + strlen(seqmod2));
          Read_Qual2(seqmod2,qual,Qualdistr2,gen);
          ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod2,qual);
          memset(qual, 0, sizeof(qual)); 
        }
      }
      else{
        //std::cout << "non flag "<< std::endl;
        Read_Qual2(seqmod,qual,Qualdistr1,gen);
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod,qual);
        memset(qual, 0, sizeof(qual));
        DNA_complement(seqmod2);
        std::reverse(seqmod2, seqmod2 + strlen(seqmod2));
        Read_Qual2(seqmod2,qual,Qualdistr2,gen);
        ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod2,qual);
        memset(qual, 0, sizeof(qual));        
      }
    }
    start_pos += fraglength;
    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
  }
  //std::cout << "out of while" << std::endl;
  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

void* Create_pe_threads(const char *fastafile,faidx_t *seq_ref,const char *output,bool Adapt_flag,const char* Adapter_1,const char* Adapter_2,int thread_no,int chr_total){
  //Loading in an creating my objects for the sequence files.
  
  int chr_idx = 0;
  int chr_no = 0;

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafq_pe_thread struct_for_threads[nthreads];

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
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].Adapter_2 = Adapter_2;
      std::cout << "chr_no" << i << std::endl;
    }
    
    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,Fafq_thread_pe_run,&struct_for_threads[j]);
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
    //srand48(start_pos+std::time(nullptr)); //you could make srand48(start_pos + seed_input)
   
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    /*std::cout << "_-------------------" << std::endl;
    std::cout << "before" << std::endl;
    std::cout << " fragment length " << fraglength << std::endl;
    std::cout << "fails" << std::endl;*/
    
    // int readlength = sizearray[SizeDist[1](gen)];
    int stop = start_pos+(int) fraglength;
    
    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+start_pos-1,150);
    }
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+start_pos-1,150);
    }
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+start_pos-1,fraglength);
    }

    //removes reads with NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += fraglength+1;} // readlength + 1
    else{
      Ill_err(seqmod,Error,gen);

      if(struct_obj->Adapter_flag==true){  
        if (fraglength < 150){ // for reads to be added with adapter

          char read1N[lib_size+1];
          char read1[lib_size + 1];

          memset(read1N,'N',lib_size);
          read1N[lib_size]='\0';

          //Copies sequence into both reads
          strcpy(read1, seqmod);
          std::strcat(read1,struct_obj->Adapter_1);

          if (strlen(read1) >= 150){
            std::strncpy(read1N, read1, 150);
          }
          else{
            //std::cout << "else " << std::endl;
            //std::cout << "read1 v2 \n" << read1 << std::endl;
            std::strncpy(read1N, read1, strlen(read1)); //copy read+adapter into N...
            //std::cout << read1N << std::endl;
          }
          
          //Creating read quality
          Read_Qual2(read1N,qual,Qualdistr1,gen);
          ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,read1N,qual);
          memset(qual, 0, sizeof(qual));    
        }
        else{
          //std::cout << "else loop" << std::endl;
          Read_Qual2(seqmod,qual,Qualdistr1,gen);
          ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod,qual);
          memset(qual, 0, sizeof(qual));
        }
      }
      else{
        //std::cout << "non flag "<< std::endl;
        Read_Qual2(seqmod,qual,Qualdistr1,gen);
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+fraglength-1,fraglength,seqmod,qual);
        memset(qual, 0, sizeof(qual));      
      }
    }
    start_pos += fraglength;
    memset(seqmod, 0, sizeof seqmod);
  }
  //std::cout << "out of while" << std::endl;
  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
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
    gz1 = gzopen("threads_test_se.fq.gz","wb");
    for(int i=0;i<nthreads;i++){
      gzwrite(gz1,struct_for_threads[i].fqresult_r1->s,struct_for_threads[i].fqresult_r1->l);
      struct_for_threads[i].fqresult_r1->l = 0;
    }
    gzclose(gz1);
  }
  else if (std::strcmp(output, "fq") == 0){
    FILE *fp1;
    fp1 = fopen("threads_test_se.fq","wb");
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

  bool Adapt_flag;
  const char* PlatformType;
  PlatformType = "SE";

  int thread_to_run = 1;
  if (std::strcmp(argv[2], "False") == 0 || std::strcmp(argv[2], "false") == 0 || std::strcmp(argv[2], "F") == 0){
    Adapt_flag = false;
    const char* Adapter_1 = NULL;
    const char* Adapter_2 = NULL;
    //Create_pe_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,Adapter_2,thread_to_run,chr_total);
    Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,thread_to_run,chr_total);
  }
  else if (std::strcmp(argv[2], "True") == 0 || std::strcmp(argv[2], "true") == 0 || std::strcmp(argv[2], "T") == 0){
    Adapt_flag = true;
    const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    const char* Adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
    //Create_pe_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,Adapter_2,thread_to_run,chr_total);
    Create_se_threads(fastafile,seq_ref,argv[1],Adapt_flag,Adapter_1,thread_to_run,chr_total);
  }
  //fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
}

// g++ SimulAncient_func.cpp FaFq_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
// ./a.out fq F
