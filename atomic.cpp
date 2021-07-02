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

/*
struct sum_runner_struct {
  int answer;
  int limit;
};

std::atomic<int> counter;

void* sum_runner(void* arg){
    struct sum_runner_struct *arg_struct = (struct sum_runner_struct*) arg;

    for(int i = 0; i <= arg_struct ->limit; i++){
        std::cout << "i " << i << std::endl;
        counter += i;
        std::cout << " count " << counter << std::endl;
    }

    //setting the answer to the sum from the calculated sum from the loop
    std::cout << "answer " << counter << std::endl;
    arg_struct -> answer = counter;

}


int main(int argc,char **argv){
    //creating an array with the arguments to create multiple threads;
    int num_args = 2;
    struct sum_runner_struct args[num_args];
    pthread_t tids[num_args];

    pthread_t mythreads[num_args];
    sum_runner_struct struct_for_threads[num_args];

    //Creating the threads 
    for (int i =0; i<=num_args; i++){
        args[i].limit = atoll(argv[i]);

        pthread_attr_t attr;

        pthread_attr_init(&attr);
        pthread_create(&tids[i],&attr,sum_runner, &args[i]);
    }

    // Wait until the threads are done before joining
    for (int i = 0;i<num_args;i++){
        pthread_join(tids[i],NULL);
    }

    std::cout << "final " << counter << std::endl;
}*/


void* full_genome_test(faidx_t *seq_ref,int chr_total){
    
  int genome_size = 0;
  int chr_sizes[chr_total];
  const char *chr_names[chr_total];
  int chr_size_cumm[chr_total];

  for (size_t i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    std::cout << "genomesize " << genome_size << std::endl;
    genome_size += chr_len;
    chr_size_cumm[i] = genome_size;
    std::cout << "genomesize 2 " << genome_size << std::endl;
  }

  std::cout << genome_size << std::endl;
  std::cout << "----------" << std::endl;
  char* genome = (char*) malloc(genome_size);
  
  for (size_t i = 0; i < chr_total; i++){

    pthread_mutex_lock(&data_mutex);
    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    pthread_mutex_unlock(&data_mutex);
    std::cout << " size " << chr_sizes[i] << std::endl;
    std::cout << " name " << chr_names[i] << std::endl;
    std::cout << " cummul " << chr_size_cumm[i] << std::endl;
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
  int chr_start;
  for (size_t i = 0; i < 10; i++){
    std::cout << "---------------------" << std::endl;
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    srand48(fraglength+std::time(nullptr));
    rand_start = lrand48() % (genome_size-fraglength-1);
    std::cout << "start " << rand_start << std::endl;
    std::cout << "frag " << fraglength << std::endl;
    
    int chr_idx = 0;
    while (rand_start > chr_size_cumm[chr_idx]){chr_idx++;}
    std::cout << "chr_idx " << chr_idx << std::endl;
    std::cout << chr_names[chr_idx] << std::endl;

    // extracting corresponding coordinates after chromosome concatenations
    if (chr_idx == 0){chr_start = rand_start;}
    else{chr_start = rand_start-chr_size_cumm[chr_idx-1];}

    fprintf(stderr,">%s:%d-%d_length:%d\n",chr_names[chr_idx],chr_start,chr_start+fraglength-1,fraglength);
    std::cout << "start " << chr_start << " end  " << chr_start+fraglength-1 << " " << chr_size_cumm[chr_idx] << std::endl;
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
