#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>//for printing time

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
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

void printTime(FILE *fp){
  time_t rawtime;
  struct tm * timeinfo; 
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}

pthread_mutex_t Fq_write_mutex;


char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],int chr_size_cumm[]){
  
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      sprintf(genome+strlen(genome),data); 
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

std::atomic<float> current_cov_atom(0.0);
std::atomic<int> current_reads_atom(0);
std::atomic<size_t> D_total(0);
std::atomic<size_t> nread_total(0);

// ---------------------- SINGLE-END ---------------------- //

struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;
  int* sizearray;
  std::discrete_distribution<> *SizeDist;
  double* Qualfreq;
  float current_cov;
  size_t current_read_no;
  int threadseed;
  float cov;
  size_t reads;
  BGZF *bgzf;
  const char* Adapter_flag;
  const char* Adapter_1;
};


void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //std::cout << "se run "<< std::endl;
  time_t t4=time(NULL);
  // creating random objects for all distributions.
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno;
  //std::random_device rd;
  //std::default_random_engine gen(loc_seed);//gen(struct_obj->threadseed); //struct_obj->seed+struct_obj->threadno
  
  size_t genome_len = strlen(struct_obj->genome);

  //coverage method2
  char seqmod[1024] = {0};
  char seqmod2[1024] = {0};
  char read[1024] = {0};
  char readadapt[1024] = {0};
  // for the coverage examples
  float cov = struct_obj -> cov;
  int reads = struct_obj -> reads;

  //float cov_current = 0;
  size_t rand_start;
  //int nread = 0;

  char qual[1024] = "";
  //int D_i = 0;
  int localread = 0;
  int iter = 0;
  int current_reads_atom = 0;
  //while (current_cov_atom < cov) {
  //while (current_reads_atom < reads){
  while (current_reads_atom < reads){
    int fraglength = 100; //(int) struct_obj->sizearray[struct_obj->SizeDist[1](gen)]; //100; //  //150; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
    //double a = ((double)rand_r(&loc_seed)/ RAND_MAX);
    //double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    //fprintf(stderr,"THREAD NO %d \t seed %lf\n",struct_obj->threadno,a);
    //srand48(1+fraglength+iter);
    //fprintf(stderr,"FRAGMENT LENGTH %d \n",fraglength);
    //fprintf(stderr,"\t-> Seed used: %d \n",seed+fraglength+iter+D_total+std::time(nullptr));

    rand_start = genome_len-100000; //rand_val * (genome_len-fraglength)-1; //lrand48() % (genome_len-fraglength-1);
    //fprintf(stderr,"RANDOM START LENGTH %d \n",rand_start);
    //identify the chromosome based on the coordinates from the cummulative size array
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}

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
    //srand48(seed+iter); 
    int rand_id = 100;//rand_val * fraglength-1;
    
    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seqmod,'N');
    pch2 = strrchr(seqmod,'N');
    //if (pch != NULL){continue;}
    if ((int )(pch-seqmod+1) == 1 && (int)(pch2-seqmod+1)  == strlen(seqmod)){continue;}
    else{
      int seqlen = strlen(seqmod);

      D_total += fraglength-1; //update the values for when calculating the coverage
      nread_total ++;
      
      //for (int j = 0; j < fraglength; j++){D_total += 1;}

      // SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097,std::time(nullptr));

      int strand = 1;//rand() % 2;

      // FASTQ FILE
      if (strand == 0){
        DNA_complement(seqmod2);
        reverseChar(seqmod2);
        //SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097);
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod2);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          strncpy(readadapt, read, 150);
          
          // Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq);
            
          ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n",struct_obj->threadno, rand_id,0,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,readadapt);
        }
        else if (struct_obj->Adapter_flag == "false"){
          
          // Read_Qual_new(seqmod2,qual,loc_seed,struct_obj->Qualfreq);
          
          ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n",struct_obj->threadno, rand_id,0,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,seqmod2);
        }
      }
      else if (strand == 1){
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          strncpy(readadapt, read, 150);
          
          Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq);      
          
          ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,1,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,readadapt,qual);
        }
        else if (struct_obj->Adapter_flag == "false"){
          
          Read_Qual_new(seqmod,qual,loc_seed,struct_obj->Qualfreq);
          
          ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,1,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,seqmod,qual);
        }
      }        

      /*if (struct_obj->fqresult_r1->l > 30000000){
        fprintf(stderr,"\t Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
        pthread_mutex_lock(&Fq_write_mutex);
        bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
        pthread_mutex_unlock(&Fq_write_mutex);
        struct_obj->fqresult_r1->l =0;
      }*/
    
      //int64_t offs[2];
      //assert(0==bgzf_flush(struct_obj->bgzf));
      //offs[0] = bgzf_tell(struct_obj->bgzf);
      memset(qual, 0, sizeof(qual));  
      //nread++;
      //fprintf(stderr,"Number of reads %d \n",nread);
    }
    
    current_cov_atom = (float) D_total / genome_len;
    current_reads_atom = (size_t) nread_total;
    struct_obj->current_cov = current_cov_atom; //why do we need this
    struct_obj->current_read_no = current_reads_atom;

    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
    chr_idx = 0;
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
    iter++;
    localread++;
    //std::cout << "currect coverage " << current_cov_atom << std::endl;
    //std::cout << "currect number of reads " << current_reads_atom << std::endl;
  }
  /*if (struct_obj->fqresult_r1->l > 0){
    fprintf(stderr,"\t last Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
    pthread_mutex_lock(&Fq_write_mutex);
    bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
    pthread_mutex_unlock(&Fq_write_mutex);
    struct_obj->fqresult_r1->l =0;
  }*/

  //delete[] struct_obj->Qualfreq;
  //delete[] struct_obj->sizearray;
  //consider freeing these after the join operator
  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  fprintf(stderr,"\t number of reads generated by this thread %d \n",localread);
  //fprintf(stderr, "\t[ALL done] walltime spend in thread %d =  %.2f sec\n", struct_obj->threadno, (float)(time(NULL) - t4));  
  //std::cout << "thread done" << std::endl;
  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, float coverage, int reads,const char* Adapt_flag,const char* Adapter_1){
  time_t t3=time(NULL);
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];
  
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);

  if (genome_data != NULL){
    fprintf(stderr,"\t-> Full genome function run!\n");
    fprintf(stderr,"\t-> Full genome size %lu \n",strlen(genome_data));
  
    
    //std::cout << " genome length " << genome_len << std::endl;
    Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

    BGZF *bgzf;
    const char* filename = "chr22_out.fq.gz";
    const char *mode = "r";
	  //bgzf = bgzf_open(filename, mode);
    bgzf = bgzf_open(filename,"wb"); 
    
    //if(nthreads>1){bgzf_mt(bgzf,nthreads,256);}
    bgzf_mt(bgzf,nthreads,256);

    std::discrete_distribution<> SizeDist[2]; 
    std::ifstream infile2("Size_dist/Size_freq.txt");
    Size_freq_dist(infile2,SizeDist,seed);//struct_obj->threadseed //creates the distribution of all the frequencies
    infile2.close();

    // Creates the random lengths array and distributions //
    std::ifstream infile("Size_dist/Size_freq.txt");
    int* sizearray = Size_select_dist(infile);
    infile.close();

    // READ QUAL ARRAY
    double* Qual_freq_array = new double[6000];
    Qual_freq_array = Qual_array(Qual_freq_array,"/home/wql443/WP1/SimulAncient/Qual_profiles/Acc_freq1.txt");
  
    //initialzie values that should be used for each thread
    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;
      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].sizearray = sizearray;
      struct_for_threads[i].SizeDist = SizeDist;
      struct_for_threads[i].Qualfreq = Qual_freq_array;
      //struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      struct_for_threads[i].cov = coverage;
      struct_for_threads[i].reads = reads;
      struct_for_threads[i].bgzf = bgzf;
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      
      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      /*for(int jj=0;jj<chr_total+1;jj++)
      struct_for_threads[i].size_cumm[jj]= chr_size_cumm[jj];*/
      
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    // free(catString);
    // delete[] sizearray;
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    // pthread_create( thread ID, struct for thread config, pointer to the function which new thread starts with, 
    // it returns a pointer to void data which is the 4th argument 
    fprintf(stderr,"Creating a bunch of threads\n"); 
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Fafq_thread_se_run,&struct_for_threads[i]);
    }
    fprintf(stderr,"Done Creating a bunch of threads\n");
    fflush(stderr);

    for (int i = 0; i < nthreads; i++)
    {  
      fprintf(stderr,"joing threads\n");fflush(stderr);
      // PRINT 
      pthread_join(mythreads[i],NULL);
      //fprintf(stderr, "\t[ANDET STED] walltime used for join =  %.2f sec\n", (float)(time(NULL) - t3));  
    }
    bgzf_close(bgzf);
    //bgzf_close(*bgzf);
    
    
    //for(int i=0;i<nthreads;i++){delete struct_for_threads[i].fqresult_r1 -> s;} //ERROR SUMMARY: 9 errors from 5 contexts (suppressed: 0 from 0 )
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      //ks_release(struct_for_threads[i].fqresult_r1);
      delete struct_for_threads[i].fqresult_r1;
      //delete[] struct_for_threads[i].Qualfreq;
      //delete[] struct_for_threads[i].sizearray;
      //free(struct_for_threads[i].fqresult_r1); // ERROR SUMMARY: 8 errors from 4 contexts (suppressed: 0 from 0)
    }
    delete[] sizearray;
    delete[] Qual_freq_array;
    free(genome_data);
  }
  return NULL;
}

// ------------------------------ //
int main(int argc,char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);
  // printTime(stderr);
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa chr1_3 chr1 hg19canon.fa  chr1_12   chr22 chr1_15 chr10_15  chr18_19  
  //const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  fprintf(stderr,"\t-> fasta load \n");
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  int Glob_seed = 1;
  int threads = 10;
  float cov = 1;
  size_t No_reads = 1e9;
  fprintf(stderr,"\t-> Seed used: %d with %d threads\n",Glob_seed,threads);
  fprintf(stderr,"\t-> Coverage used for simulation: %f\n",cov);
  fprintf(stderr,"\t-> Number of simulated reads: %zd\n",No_reads);
  //char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  const char* Adapt_flag = "false";
  const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
  //fprintf(stderr,"Creating a bunch of threads\n");
  Create_se_threads(seq_ref,threads,Glob_seed,cov,No_reads,Adapt_flag,Adapter_1);
  //fprintf(stderr,"Done creating a bunch of threads\n");
  //fflush(stderr);
  // free the calloc memory from fai_read
  //free(seq_ref);
  fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
}

//SimBriggsModel(seqmod, frag, L, 0.024, 0.36, 0.68, 0.0097);
// g++ SimulAncient_func.cpp atomic_fq.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
//cat test.fq | grep 'T0' | grep 'chr20' | wc -l
//valgrind --tool=memcheck --leak-check=full ./a.out

//awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' chr22_out.fq
// { time ./a.out; } 2> test.txt
// for i in {1..10}; do (time ./a.out) &>> time/cpu_2022.txt ; done
// for i in {1..10}; do (time ./a.out) &>> time/hg19c1t1.txt ; done
//for i in {1..2}; do (/usr/bin/time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U" ./hg19c1t1) &>> time/hg19c1t1.txt ; done
//for i in {1..10}; do (/usr/bin/time -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U" ./hg19c1t16) &>> time/hg19gzt2.txt ; done