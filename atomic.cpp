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
#include <mutex>        
#include <atomic>
#include <vector>

#include "SimulAncient_func.h"

pthread_mutex_t data_mutex;


char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],int chr_size_cumm[]){
  
  int genome_size = 0;
  chr_size_cumm[0] = 0;
  for (size_t i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(genome_size);
  
  for (size_t i = 0; i < chr_total; i++){

    pthread_mutex_lock(&data_mutex);
    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    pthread_mutex_unlock(&data_mutex);
    strcat(genome,data);
  }
  return genome;
}

std::atomic<float> current_cov_atom(0.0);
std::atomic<int> size_data(0);
std::atomic<int> D_total(0);

struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  char *genome; // The actual concatenated genome
  int chr_no;
  int *size;
  int *size_cumm;
  const char **names;
  const char* Ill_err;
  const char* read_err_1;
  float current_cov;
  int cov_size;
  int depth;
};

void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //std::cout << "se run "<< std::endl;

  int chr_total = struct_obj->chr_no;

  // Load in the error profiles and size distributions
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
  // -------------------------- // 
  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  int genome_len = strlen(struct_obj->genome);
  //std::cout << genome_len << std::endl;

  int start_pos = 1;
  int end_pos = genome_len; //30001000

  //coverage method2
  int* chrlencov = (int *) calloc(genome_len,sizeof(int));

  char seqmod[1024] = {0};
  char seqmod2[1024] = {0};
  // for the coverage examples
  float cov = 2;
  //float cov_current = 0;
  int rand_start;
  int fraglength;
  int nread = 0;

  char qual[1024] = "";
  //int D_i = 0;
  
  
  while (current_cov_atom < cov) {
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
    srand48(D_total+fraglength+std::time(nullptr));
    rand_start = lrand48() % (genome_len-fraglength-1);
    /*std::cout << "start " << rand_start << std::endl;
    std::cout << "frag " << fraglength << std::endl;*/

    //identify the chromosome based on the coordinates from the cummulative size array
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}
    //std::cout << "chr_idx " << chr_idx << std::endl;
    //std::cout << struct_obj->names[chr_idx] << std::endl;

    // case 1
    if (fraglength > 150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,struct_obj->genome+rand_start-1,150);
    }
    // case 2
    else if (fraglength <= 150)
    {
      strncpy(seqmod,struct_obj->genome+rand_start-1,fraglength);
    }

    //removes reads with NNN
    int rand_id = (lrand48() % (genome_len-fraglength-1))%fraglength;
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){continue;}
    else{
      for (int j = rand_start; j < rand_start+fraglength; j++){
        chrlencov[j] += 1; //adds 1 for regions, which might have reads already
        D_total += 1; // to find the total depth
        if (chrlencov[j] == 1){size_data++;} // if the value is different from 1 (more reads) then we still count that as one chromosome site with data
      }
      SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097);
      Ill_err(seqmod2,Error,gen);
      Read_Qual2(seqmod2,qual,Qualdistr1,gen);
      ksprintf(struct_obj->fqresult_r1,"@READID_%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",rand_id,
        struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
        fraglength,seqmod2,qual);
      memset(qual, 0, sizeof(qual));  
      nread++;
      //fprintf(stderr,"Number of reads %d , d_total %d, size_data %d\n",nread,D_total,size_data);
    }
    current_cov_atom = (float) D_total / (genome_len-size_data);
    //std::cout << "----------------------------" << std::endl;
    //std::cout << "current " << current_cov_atom << std::endl;
    struct_obj->current_cov = current_cov_atom; //why do we need this
    struct_obj->cov_size = size_data;
    struct_obj->depth = D_total;

    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
    chr_idx = 0;
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
  }
  //std::cout << "thread done" << std::endl;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no){
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];
  
  int chr_total = faidx_nseq(seq_ref);

  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);

  Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

  //initialzie values that should be used for each thread
  for (int i = 0; i < nthreads; i++){
    struct_for_threads[i].fqresult_r1 =new kstring_t;
    struct_for_threads[i].fqresult_r1 -> l = 0;
    struct_for_threads[i].fqresult_r1 -> m = 0;
    struct_for_threads[i].fqresult_r1 -> s = NULL;
    struct_for_threads[i].genome = genome_data;
    struct_for_threads[i].chr_no = chr_total;
    struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
    struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
    
    //declaring the size of the different arrays
    struct_for_threads[i].size = (int*)malloc(sizeof(int) * struct_for_threads[i].chr_no); 
    memcpy(struct_for_threads[i].size, chr_sizes, sizeof(chr_sizes));
    struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * struct_for_threads[i].chr_no);
    memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
    struct_for_threads[i].names = (const char**)malloc(sizeof(char) * struct_for_threads[i].chr_no);
    memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
  }

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  // pthread_create( thread ID, struct for thread config, pointer to the function which new thread starts with, 
  // it returns a pointer to void data which is the 4th argument  
  for (int i = 0; i < nthreads; i++){
    pthread_create(&mythreads[i],&attr,Fafq_thread_se_run,&struct_for_threads[i]);
  }
  
  for (int i = 0; i < nthreads; i++)
  {
    pthread_join(mythreads[i],NULL);
  }
  
  FILE *fp1;
  fp1 = fopen("lol.fq","wb");
  for(int i=0;i<nthreads;i++){
    if (struct_for_threads[i].fqresult_r1->l > 10000){
      fwrite(struct_for_threads[i].fqresult_r1->s,sizeof(char),struct_for_threads[i].fqresult_r1->l,fp1);struct_for_threads[i].fqresult_r1->l =0;
    }
  }
  /*
  for(int i=0;i<nthreads;i++){
    fprintf(fp1,"%s",struct_for_threads[i].fqresult_r1->s);
  }*/
  fclose(fp1);
  return NULL;
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
  int threads = 5;
  fprintf(stderr,"\t-> %d threads\n",threads);
  //full_genome_create(seq_ref,chr_total);
  Create_se_threads(seq_ref,threads);
}

//SimBriggsModel(seqmod, frag, L, 0.024, 0.36, 0.68, 0.0097);
// g++ SimulAncient_func.cpp atomic.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl

