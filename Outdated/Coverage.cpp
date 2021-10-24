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

struct Parsarg_for_Fafa_thread{
  int chr_idx;
  kstring_t *fqresult;
  faidx_t *seq_ref;
  //const char* Output_type;
};

//-------------- FUNCTIONS ----------------//

void coverage(float cov){
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr1.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;
  const char *chr_name = faidx_iseq(seq_ref,chr_no);
  int chr_len = 1000000;//faidx_seq_len(seq_ref,chr_name); //

  int* chrlencov = (int *) calloc(chr_len,sizeof(int));
  std::cout << "chr_len " << faidx_seq_len(seq_ref,chr_name)  << std::endl;
  char seqmod[1024] = {0};

  // Creates the random lengths array and distributions //
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 

  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  float cov_current = 0;
  int rand_start;
  int fraglength;
  int nread = 0;

  while (cov_current < cov) {
    int D_i = 0;
    int N = chr_len; // N number of sites with data > 0
    std::cout << "---------------------" << std::endl;
    std::cout << "COUNT " << N << std::endl;
    
    fraglength = (int) sizearray[SizeDist[1](gen)];
    //std::cout << "fraglenth" << fraglength << std::endl;
    srand48(fraglength+std::time(nullptr));
    rand_start = lrand48() % (chr_len-fraglength-1);
    
    /*std::cout << rand_start << std::endl;
    fprintf(stderr,"start %d length %d",rand_start,fraglength);
    std::cout << "---------" << std::endl;*/

    for (int j = rand_start; j < rand_start+fraglength; j++){chrlencov[j] += 1;}
    //for (size_t i = 10; i < 1000; i++){std::cout << chrlencov[i] << ' ';}
    std::cout << std::endl;
    for(int i = 0; i<N ; i++){
      if (chrlencov[i] == 0){N = N-1;} //trækker -1 fra da 0 jo ikke er data -> så alle positioner med værdi forskellig for 0
      else{D_i+=chrlencov[i];} // opdaterer summen D_i = D_i+j for length
    }
    nread++;
    cov_current = (float) D_i / N; // update the mean depth
    fprintf(stderr,"Number of reads %d , sum %d, count %d, coverage %f\n",nread,D_i,N,cov_current);
    fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
    //std::cout << "number reads " << nread << " current coverage "<< cov_current << std::endl;
  }
}

void random_seq(faidx_t *seq_ref,int nreads){
  // choose a random sequence -> still ned to change it so it saves the output to a single file.
  int readlength=35;
  // add size dist etc
  for(int i=0;i<nreads;i++){
    char buf[96];//assume maxlength for readid is 96bytes
    int whichref = lrand48() % faidx_nseq(seq_ref);
    fprintf(stderr,"\t-> Whichref: %d\n",whichref);
    const char *name = faidx_iseq(seq_ref,whichref);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"\t-> name: \'%s\' name_len: %d\n",name,name_len);

    int start = lrand48() % name_len;
    int stop = start+readlength;
    if(stop>name_len)
      stop = name_len;
    snprintf(buf,96,"%s:%d-%d",name,start,stop);
    fprintf(stderr,"buf: %s\n",buf);
    fprintf(stdout,"%s\n+\n",buf);
    char *data = fai_fetch(seq_ref,name,&name_len);

    for(int i=start;i<stop;i++)
      fprintf(stdout,"%c",data[i]);
    fprintf(stdout,"\n");

    for(int i=start;i<stop;i++)
      fprintf(stdout,"F");
    fprintf(stdout,"\n");
  }
}

// --------------------- //
int main(int argc,char **argv){
  coverage(4);
}

// g++ SimulAncient_func.cpp FaFa_thread_cov.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// ./a.out fa 5