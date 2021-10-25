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

    pthread_mutex_lock(&data_mutex);
    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    pthread_mutex_unlock(&data_mutex);
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
std::atomic<int> size_data(0);
std::atomic<int> D_total(0);
std::atomic<int> nread_total(0);

// ---------------------- SINGLE-END ---------------------- //

struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;
  const char* Ill_err;
  float current_cov;
  int cov_size;
  int threadseed;
  float cov;
  FILE *fp1;
  const char* Adapter_flag;
  const char* Adapter_1;

  // VCF
  bcf1_t *vcf_record;
  htsFile *in_vcf;
  bcf_hdr_t *vcf_hdr;
};


void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //std::cout << "se run "<< std::endl;
  
  // creating random objects for all distributions.
  int seed = struct_obj->threadseed+struct_obj->threadno;
  std::random_device rd;
  std::default_random_engine gen(seed);//gen(struct_obj->threadseed); //struct_obj->seed+struct_obj->threadno
  
  // -------------------------- // 
// Load in the error profiles and size distributions
  
  double** Error_2darray = create2DArray(struct_obj->Ill_err,4,280);
  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);
  
  //Free each sub-array
  for(int i = 0; i < 280; ++i) {delete[] Error_2darray[i];}
  //Free the array of pointers
  delete[] Error_2darray;
  
  // -------------------------- // 
  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_dist/Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_dist/Size_freq.txt");
  Size_freq_dist(infile2,SizeDist,seed);//struct_obj->threadseed //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //
  size_t genome_len = strlen(struct_obj->genome);
  //coverage method2
  int* chrlencov = (int *) calloc(genome_len,sizeof(int));
  char seqmod[1024] = {0};
  char seqmod2[1024] = {0};
  char read[1024] = {0};
  char readadapt[1024] = {0};
  // for the coverage examples
  float cov = struct_obj -> cov;
  //float cov_current = 0;
  size_t rand_start;
  size_t rand_start_tmp;
  int nread = 0;

  /*srand48(seed+std::time(nullptr));
  rand_start_tmp = lrand48() % (genome_len-1) + seed;
  fprintf(stderr,"temp start %d\n",rand_start_tmp);*/

  const char *chrID = bcf_hdr_id2name(struct_obj->vcf_hdr, struct_obj->vcf_record->rid);

  hts_idx_t *idx = bcf_index_load("/home/wql443/WP1/Alba/HLA_test.bcf.csi");
  if(!idx){printf("Null index\n");}

  //int D_i = 0;
  /*char region[128];
    snprintf(region, sizeof(region), "%s:%d-%d", chrID,1000000,2000000);
    //fprintf(stderr,"Chromosomal region from vcf -> %s \t chromosome: %s\n",region,struct_obj->names[0]);
    hts_itr_t *itr = bcf_itr_querys(idx, struct_obj->vcf_hdr, region);
    if(!itr){printf("Null iterator\n");}
    //fprintf(stderr,"-------------- \n");
    while ((bcf_itr_next(struct_obj->in_vcf, itr, struct_obj->vcf_record)) == 0) {
      bcf_unpack((bcf1_t*)struct_obj->vcf_record, BCF_UN_ALL);
      if (struct_obj->vcf_record->n_allele > 1){
        printf("id: %d, pos: %d\n", struct_obj->vcf_record->rid, struct_obj->vcf_record->pos);
        //printf("id: %d, pos: %d, len: %d, allele: %d, allele: %s\n", struct_obj->vcf_record->rid, struct_obj->vcf_record->pos, struct_obj->vcf_record->rlen,struct_obj->vcf_record->n_allele,struct_obj->vcf_record->d.allele[0]);
      }
    }
    bcf_itr_destroy(itr);
    */
  int iter = 0;
  while (current_cov_atom < cov) {
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //150; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
    char region[128];
    snprintf(region, sizeof(region), "%s:%d-%d", chrID,rand_start-1,rand_start-1+fraglength);
    //fprintf(stderr,"Chromosomal region from vcf -> %s \t chromosome: %s\n",region,struct_obj->names[0]);
    hts_itr_t *itr = bcf_itr_querys(idx, struct_obj->vcf_hdr, region);
    if(!itr){printf("Null iterator\n");}
    //fprintf(stderr,"-------------- \n");
    while ((bcf_itr_next(struct_obj->in_vcf, itr, struct_obj->vcf_record)) == 0) {
      bcf_unpack((bcf1_t*)struct_obj->vcf_record, BCF_UN_ALL);
      //if (struct_obj->vcf_record->n_allele > 1){
        //printf("id: %d, pos: %d\n", struct_obj->vcf_record->rid, struct_obj->vcf_record->pos);
        //printf("id: %d, pos: %d, len: %d, allele: %d, allele: %s\n", struct_obj->vcf_record->rid, struct_obj->vcf_record->pos, struct_obj->vcf_record->rlen,struct_obj->vcf_record->n_allele,struct_obj->vcf_record->d.allele[0]);
      //}
    }
    bcf_itr_destroy(itr);
    
    srand48(seed+fraglength+iter);
    //srand48(seed+fraglength+iter+D_total+std::time(nullptr)); //D_total+fraglength //+std::time(nullptr)
    //fprintf(stderr,"\t-> Seed used: %d \n",seed+fraglength+iter+D_total+std::time(nullptr));

    rand_start = lrand48() % (genome_len-fraglength-1) + seed;

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
    srand48(seed+iter); 
    int rand_id = (lrand48() % (genome_len-fraglength-1))%fraglength;
    
    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seqmod,'N');
    pch2 = strrchr(seqmod,'N');
    //if (pch != NULL){continue;}
    if ((int )(pch-seqmod+1) == 1 && (int)(pch2-seqmod+1)  == strlen(seqmod)){continue;}
    else{
      for (size_t j = rand_start; j < rand_start+fraglength; j++){
        //std::cout<< "j" << j << std::endl;
        chrlencov[j] += 1; //adds 1 for regions, which might have reads already
        
        D_total += 1; // to find the total depth
        if (chrlencov[j] == 1){size_data++;} // if the value is different from 1 (more reads) then we still count that as one chromosome site with data
      }
      
      SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097,std::time(nullptr));

      int strand = rand() % 2;

      // FASTQ FILE
      if (strand == 0){
        DNA_complement(seqmod2);
        reverseChar(seqmod2);
        //SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097);
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod2);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          Ill_err(read,Error,gen);
          strncpy(readadapt, read, 150);          
          ksprintf(struct_obj->fqresult_r1,">T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,readadapt);
        }
        else if (struct_obj->Adapter_flag == "false"){
          Ill_err(seqmod2,Error,gen);
          ksprintf(struct_obj->fqresult_r1,">T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,seqmod2);
        }
      }
      else if (strand == 1){
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod2);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          Ill_err(read,Error,gen);
          strncpy(readadapt, read, 150);
          ksprintf(struct_obj->fqresult_r1,">T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,readadapt);
        }
        else if (struct_obj->Adapter_flag == "false"){
          Ill_err(seqmod2,Error,gen);
          ksprintf(struct_obj->fqresult_r1,">T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength,seqmod2);
        }
      }        
      //struct_obj->fqresult_r1->l =0;
      if (struct_obj->fqresult_r1->l > 100000){
        pthread_mutex_lock(&Fq_write_mutex);
        fwrite(struct_obj->fqresult_r1->s,sizeof(char),struct_obj->fqresult_r1->l,struct_obj->fp1);
        pthread_mutex_unlock(&Fq_write_mutex);
        struct_obj->fqresult_r1->l =0;
      }

      nread++;
      //fprintf(stderr,"Number of reads %d \n",nread);
    }
    current_cov_atom = (float) D_total / genome_len;

    struct_obj->current_cov = current_cov_atom; //why do we need this
    struct_obj->cov_size = size_data;

    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
    chr_idx = 0;
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
    iter++;
    //std::cout << current_cov_atom << "\t" << cov << std::endl;
  }
  hts_idx_destroy(idx);
  //consider freeing these after the join operator
  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  free(chrlencov);
  
  //std::cout << "thread done" << std::endl;
  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, float coverage,const char* Adapt_flag,const char* Adapter_1){
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];
  
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);

      // -------- VCF INFO -------- //
  bcf1_t *vcf_record = bcf_init();
  htsFile *in = bcf_open("/home/wql443/WP1/Alba/HLA_test.bcf", "r");
  if (!in){printf("Null file pointer\n");}
  bcf_hdr_t *hdr = bcf_hdr_read(in);

  if (genome_data != NULL){
    fprintf(stderr,"\t-> Full genome function run!\n");
    fprintf(stderr,"\t-> Full genome size %lu \n",strlen(genome_data));

  //std::cout << " genome length " << genome_len << std::endl;
    Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

    /*const char *one = "chr22_out.";
    char* catString = (char*)malloc(sizeof(char) * strlen(one) + 1); //malloc(strlen(one)+strlen(two)+1);
    strcpy(catString, one);
    strcat(catString, output);*/
    //DOING THIS FUCKED UP WITH MORE VALGRIND ERRORS
    
    FILE *fp1;
    fp1 = fopen("chr6_out.fa","wb");

    // Error profiles
    double** Error_2darray = create2DArray("/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt",4,280);
    std::discrete_distribution<> Error[280];
    Seq_err(Error_2darray,Error,280);
    
    //Free each sub-array
    for(int i = 0; i < 280; ++i) {delete[] Error_2darray[i];}
    //Free the array of pointers
    delete[] Error_2darray;

    std::discrete_distribution<> SizeDist[2]; 
    std::ifstream infile2("Size_dist/Size_freq.txt");
    Size_freq_dist(infile2,SizeDist,seed);//struct_obj->threadseed //creates the distribution of all the frequencies
    infile2.close();

    // Creates the random lengths array and distributions //
    std::ifstream infile("Size_dist/Size_freq.txt");
    int* sizearray = Size_select_dist(infile);
    infile.close();

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
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].cov = coverage;
      struct_for_threads[i].fp1 = fp1;
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;

      struct_for_threads[i].vcf_record = vcf_record;
      struct_for_threads[i].in_vcf = in;
      struct_for_threads[i].vcf_hdr = hdr;
      
      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    //free(catString);
    delete[] sizearray;
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
    fclose(fp1);
    
    
    //for(int i=0;i<nthreads;i++){delete struct_for_threads[i].fqresult_r1 -> s;} //ERROR SUMMARY: 9 errors from 5 contexts (suppressed: 0 from 0 )
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      //ks_release(struct_for_threads[i].fqresult_r1);
      delete struct_for_threads[i].fqresult_r1;
      //free(struct_for_threads[i].fqresult_r1); // ERROR SUMMARY: 8 errors from 4 contexts (suppressed: 0 from 0)
    }
 
    free(genome_data);

    //hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    bcf_destroy1(vcf_record);
    bcf_close(in);
  }
  return NULL;
}

// ------------------------------ //
int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa  hg19canon.fa  chr1_12   chr22 chr1_15 chr10_15
  //const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr6.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  fprintf(stderr,"\t-> fasta load \n");
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  
  int seed = 1;
  int threads = 1;
  float cov = 1;
  fprintf(stderr,"\t-> Seed used: %d with %d threads\n",seed,threads);

  //char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  const char* Adapt_flag = "false";
  const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";

  Create_se_threads(seq_ref,threads,seed,cov,Adapt_flag,Adapter_1);

  // free the calloc memory from fai_read
  //free(seq_ref);
  fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks

}

//SimBriggsModel(seqmod, frag, L, 0.024, 0.36, 0.68, 0.0097);
// g++ SimulAncient_func.cpp atomic_VCF2.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -o vcf
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
//cat test.fq | grep 'T0' | grep 'chr20' | wc -l
//valgrind --tool=memcheck --leak-check=full ./a.out

//awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' chr22_out.fq
// grep 
