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

//#include "SimulAncient_func.h"

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
  const char* Adapter_1;
  const char* Adapter_2;
};

//-------------- FUNCTIONS ----------------//
double** create2DArray(const char* filename,int width,int height){
  /* create 2d objects for a given error profile, with rows being position in read for specific nt
  and the cells being the frequency values */

  std::ifstream infile(filename);
  double** array2D = 0;
  array2D = new double*[height];
  for (int h = 0; h < height; h++){
    array2D[h] = new double[width];
    for (int w = 0; w < width; w++){
      infile >> array2D[h][w];}
  }
  infile.close();
  return array2D;
}

void Read_Qual2(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  /* creating the nucleotide quality string for fastq format and bam format using char array*/

  const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};

  int read_length = strlen(seq);

  //the line offset for the distribution *Dist created from the 600*8 2Darray
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  // Dist[row_idx](gen) returns values between 0-7 (columns) sampling one out of 8 nt qual from the 
  // for each line (read positions) in the read using the error profile created using the 2d array 
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        strncat(qual, nt_qual[Dist[row_idx](gen)], 1);
        break;
      case 'T':
      case 't':
        strncat(qual, nt_qual[Dist[row_idx + Tstart](gen)], 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, nt_qual[Dist[row_idx + Gstart](gen)], 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, nt_qual[Dist[row_idx + Cstart](gen)], 1);
        break;
      case 'N':
        strncat(qual, nt_qual[0], 1);;
        break;
    }
  }
}

void DNA_complement(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'T';
        break;
      case 'G':
      case 'g':
        *seq = 'C';
        break;
      case 'C':
      case 'c':
        *seq = 'G';
        break;
      case 'T':
      case 't':
        *seq = 'A';
        break;  
    }
    ++seq;
  }
}

void Qual_dist(double** Array2d,std::discrete_distribution<> dist[],int size){
  /* creating a discrete distribution of nucleotide qualitie for each line in 2D array */
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }
}

void Seq_err(double** Array2d,std::discrete_distribution<> nt_sub[],int size){
  /* Similar to qual_dist creating a nucleotide distribution for each line in 2D array*/
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],
                                    Array2d[row_idx][1],
                                    Array2d[row_idx][2],
                                    Array2d[row_idx][3]});  
    nt_sub[row_idx] = d;
  }
}

void Ill_err(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  /* Similar to Read_Qual2 but for creating substitutions in the string */
  int read_length = strlen(seq);
  int Tstart = 70;
  int Gstart = 140;
  int Cstart = 210;
  const char LookUp_nt[4] = {'A','T','G','C'};
  
  for (int nt_idx = 0; nt_idx < read_length; nt_idx++){
    switch(seq[nt_idx]){
      case 'A':
      case 'a':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx](gen)];
        break;
      case 'T':
      case 't':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Tstart](gen)];
        break;  
      case 'G':
      case 'g':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Gstart](gen)];
        break;
      case 'C':
      case 'c':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Cstart](gen)];
        break;
    }
  }
}

// --------------------- //
void* Fafq_thread_run(void *arg){
  Parsarg_for_Fafq_thread *struct_obj = (Parsarg_for_Fafq_thread*) arg;

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
  char qual[1024] = "";
  
  while(start_pos <= end_pos){
    //std::cout << "while loop" << std::endl;
    int readlength = drand48()*(70.0-30.0)+30.0; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    int stop = start_pos+(int) readlength;
    
    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);
    //std::cout << "sequence\t" << seqmod << std::endl;
    
    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += readlength + 1;}
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
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength,readlength,read1N,qual);
        memset(qual, 0, sizeof(qual));    

        Read_Qual2(read2N,qual,Qualdistr2,gen);
        ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength,readlength,read2N,qual);
        memset(qual, 0, sizeof(qual));    
      }
      else{
        //std::cout << "non flag "<< std::endl;
        Read_Qual2(seqmod,qual,Qualdistr1,gen);
        ksprintf(struct_obj->fqresult_r1,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength,readlength,seqmod,qual);
        memset(qual, 0, sizeof(qual));
        DNA_complement(seqmod);
        std::reverse(seqmod, seqmod + strlen(seqmod));
        Read_Qual2(seqmod,qual,Qualdistr2,gen);
        ksprintf(struct_obj->fqresult_r2,"@%s:%d-%d_length:%d\n%s\n+\n%s\n",chr_name,start_pos,start_pos+readlength,readlength,seqmod,qual);
        memset(qual, 0, sizeof(qual));        
      }
    }
    start_pos += readlength + 1;
    memset(seqmod, 0, sizeof seqmod);
    //std::cout << "start pos "<< start_pos << std::endl;
    //std::cout << "end pos "<< end_pos << std::endl;
  }
  //std::cout << "out of while" << std::endl;
  //for(int i=0;i<100000;i++){ksprintf(struct_obj->fqresult,"@%s_%i_%i\nCGTGA\n+\nIIIII\n",chr_name,chr_len,idx);}
  //now we have generated 1mio fq reads.
  std::cout << chr_name << " chromosome done" << std::endl;
  return NULL;
}

int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  //fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));

  bool Adapt_flag;
  const char* Adapter_1;
  const char* Adapter_2;
  if (std::strcmp(argv[2], "False") == 0 || std::strcmp(argv[2], "false") == 0 || std::strcmp(argv[2], "F") == 0){
    Adapt_flag = false;
  }
  else if (std::strcmp(argv[2], "True") == 0 || std::strcmp(argv[2], "true") == 0 || std::strcmp(argv[2], "T") == 0){
    Adapt_flag = true;
    Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    Adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
  }

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=chr_no;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fafq_thread struct_for_threads[nthreads];

  //initialzie values that should be used for each thread
  for(int i=0;i<nthreads;i++){
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
  }

  //launch all worker threads
  for(int i=0;i<nthreads;i++){
    pthread_attr_t attr;
    pthread_attr_init(&attr);

    pthread_create(&mythreads[i],&attr,Fafq_thread_run,&struct_for_threads[i]);
  }
    
  //right now 22 thread are running at the same time
  //now join, this means waiting for all worker threads to finish
  for(int i=0;i<nthreads;i++){
    pthread_join(mythreads[i],NULL);
  }

  pthread_mutex_destroy(&data_mutex);    
    //now all are done and we only write in main program.
    //now write everything

  if (std::strcmp(argv[1], "gz") == 0){
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
  else if (std::strcmp(argv[1], "fq") == 0){
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
}

// g++ TK_threads.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// cat test.fa | grep '>' | cut -c 1-6 | sort -u
