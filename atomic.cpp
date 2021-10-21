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

std::atomic<int> C_total_1(0);
std::atomic<int> C_to_T_1(0);
std::atomic<int> C_total_2(0);
std::atomic<int> C_to_T_2(0);
std::atomic<int> C_total_3(0);
std::atomic<int> C_to_T_3(0);

std::atomic<int> G_total_1(0);
std::atomic<int> G_to_A_1(0);
std::atomic<int> G_total_2(0);
std::atomic<int> G_to_A_2(0);
std::atomic<int> G_total_3(0);
std::atomic<int> G_to_A_3(0);

// ---------------------- SINGLE-END ---------------------- //

struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  const char* output;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;
  const char* Ill_err;
  const char* read_err_1;
  float current_cov;
  int cov_size;
  int depth;
  int threadseed;
  float cov;
  FILE *fp1;
  gzFile gz1;
  const char* Adapter_flag;
  const char* Adapter_1;
};


void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //std::cout << "se run "<< std::endl;
  
  // creating random objects for all distributions.
  int seed = struct_obj->threadseed+struct_obj->threadno;
  std::random_device rd;
  std::default_random_engine gen(seed);//gen(struct_obj->threadseed); //struct_obj->seed+struct_obj->threadno

  // Load in the error profiles and size distributions
  std::ifstream file(struct_obj->read_err_1);
  int Line_no = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();
  size_t lib_size = Line_no/4; 

  // loading in error profiles for nt substitions and read qual creating 2D arrays
  double** R1_2Darray = create2DArray(struct_obj->read_err_1,8,Line_no);
  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);

  //Free each sub-array
  for(int i = 0; i < Line_no; ++i) {delete[] R1_2Darray[i];}
  //Free the array of pointers
  delete[] R1_2Darray;
  
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
  int nread = 0;

  char qual[1024] = "";
  //int D_i = 0;
  
  int iter = 0;
  while (current_cov_atom < cov) {
    int fraglength = (int) sizearray[SizeDist[1](gen)]; //150; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
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
        //std::cout << D_total << std::endl;
        if (chrlencov[j] == 1){size_data++;} // if the value is different from 1 (more reads) then we still count that as one chromosome site with data
      }
      
      int seqlen = strlen(seqmod);
      if (seqmod[0]=='C' || seqmod[0]=='c'){C_total_1++;}
      if (seqmod[1]=='C' || seqmod[1]=='c'){C_total_2++;}
      if (seqmod[2]=='C' || seqmod[2]=='c'){C_total_3++;}

      if (seqmod[seqlen-1]=='G' || seqmod[seqlen-1]=='g'){G_total_1++;}
      if (seqmod[seqlen-2]=='G' || seqmod[seqlen-2]=='g'){G_total_2++;}
      if (seqmod[seqlen-3]=='G' || seqmod[seqlen-3]=='g'){G_total_3++;}

      SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097,std::time(nullptr));
      
      if ((seqmod2[0] == 'T' || seqmod2[0] == 't') && (seqmod[0] == 'C' || seqmod[0] == 'c')){C_to_T_1++;}
      if ((seqmod2[1] == 'T' || seqmod2[1] == 't')  && (seqmod[1] == 'C' || seqmod[1] == 'c')){C_to_T_2++;}
      if ((seqmod2[2] == 'T' || seqmod2[2] == 't')  && (seqmod[2] == 'C' || seqmod[2] == 'c')){C_to_T_3++;}
        
      if ((seqmod2[seqlen-1] == 'A' || seqmod2[seqlen-1] == 'a') && (seqmod[seqlen-1] == 'G' || seqmod[seqlen-1] == 'g')){G_to_A_1++;}
      if ((seqmod2[seqlen-2] == 'A' || seqmod2[seqlen-2] == 'a')  && (seqmod[seqlen-2] == 'G' || seqmod[seqlen-2] == 'g')){G_to_A_2++;}
      if ((seqmod2[seqlen-3] == 'A' || seqmod2[seqlen-3] == 'a')  && (seqmod[seqlen-3] == 'G' || seqmod[seqlen-3] == 'g')){G_to_A_3++;}

      int strand = rand() % 2;

      // FASTQ FILE
      if (struct_obj->output == "fq" || struct_obj->output == "fastq" || struct_obj->output == "fq.gz" || struct_obj->output == "fastq.gz"){
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
            Read_Qual2(readadapt,qual,Qualdistr1,gen);
            
            ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,seed,
            struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
            fraglength,readadapt,qual);
          }
          else if (struct_obj->Adapter_flag == "false"){
            Ill_err(seqmod2,Error,gen);
            Read_Qual2(seqmod2,qual,Qualdistr1,gen);
            ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,seed,
            struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
            fraglength,seqmod2,qual);
          }
        }
        else if (strand == 1){
          if (struct_obj->Adapter_flag == "true"){
            strcpy(read, seqmod2);
            strcat(read,struct_obj->Adapter_1);
            //std::cout << "read " << read << std::endl;
            Ill_err(read,Error,gen);
            strncpy(readadapt, read, 150);
            Read_Qual2(readadapt,qual,Qualdistr1,gen);
            ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,seed,
            struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
            fraglength,readadapt,qual);
          }
          else if (struct_obj->Adapter_flag == "false"){
            Ill_err(seqmod2,Error,gen);
            Read_Qual2(seqmod2,qual,Qualdistr1,gen);
            ksprintf(struct_obj->fqresult_r1,"@T%d_RID%d_S%d_%s:%d-%d_length:%d\n%s\n+\n%s\n",struct_obj->threadno, rand_id,seed,
            struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
            fraglength,seqmod2,qual);
          }
        }        
      }

      // FASTA FILE
      else if (struct_obj->output == "fa" || struct_obj->output == "fasta" || struct_obj->output == "fa.gz" || struct_obj->output == "fasta.gz"){
        //fprintf(stderr,"\t->strand\'%d\n'",strand);
        if (strand == 0){ // Reverse strand
          DNA_complement(seqmod2);
          reverseChar(seqmod2);
          //Ill_err(seqmod2,Error,gen);
          //fprintf(stderr,"\t->sequence\'%s\n'",seqmod2);
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
      }

      //struct_obj->fqresult_r1->l =0;
      if (struct_obj->fqresult_r1->l > 1000000){
        if (struct_obj->output == "fq" || struct_obj->output == "fastq" || struct_obj->output == "fa" || struct_obj->output == "fasta"){
          pthread_mutex_lock(&Fq_write_mutex);
          fwrite(struct_obj->fqresult_r1->s,sizeof(char),struct_obj->fqresult_r1->l,struct_obj->fp1);
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->fqresult_r1->l =0;
        }
        else if (struct_obj->output == "fq.gz" || struct_obj->output == "fastq.gz"){
          pthread_mutex_lock(&Fq_write_mutex);
          gzwrite(struct_obj->gz1,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->fqresult_r1->l =0;
        }
      }

      memset(qual, 0, sizeof(qual));  
      nread++;
      //fprintf(stderr,"Number of reads %d , d_total %d, size_data %d\n",nread,D_total,size_data);
    }
    current_cov_atom = (float) D_total / genome_len;
    //std::cout << "----------------------------" << std::endl;
    //std::cout << "current " << current_cov_atom << std::endl;
    struct_obj->current_cov = current_cov_atom; //why do we need this
    struct_obj->cov_size = size_data;
    struct_obj->depth = D_total;

    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);
    chr_idx = 0;
    //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
    iter++;
    //std::cout << current_cov_atom << "\t" << cov << std::endl;
  }

  //consider freeing these after the join operator
  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  free(chrlencov);
  delete[] sizearray;
  
  //std::cout << "thread done" << std::endl;
  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, float coverage, const char* output,const char* Adapt_flag,const char* Adapter_1){
  
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

    /*const char *one = "chr22_out.";
    char* catString = (char*)malloc(sizeof(char) * strlen(one) + 1); //malloc(strlen(one)+strlen(two)+1);
    strcpy(catString, one);
    strcat(catString, output);*/
    //DOING THIS FUCKED UP WITH MORE VALGRIND ERRORS
    
    FILE *fp1;
    gzFile gz1;
    if (output == "fq" || output == "fastq"){
      //FILE *fp1;  error: 'fp1' was not declared in this scope struct_for_threads[i].fp1 = fp1;
      fp1 = fopen("chr22_out.fq","wb");
    }
    else if (output == "fa" || output == "fasta"){
      //FILE *fp1;  error: 'fp1' was not declared in this scope struct_for_threads[i].fp1 = fp1;
      fp1 = fopen("chr22_out.fa","wb");
    }
    else if (output == "fq.gz" || output == "fastq.gz"){
      //FILE *fp1;  error: 'fp1' was not declared in this scope struct_for_threads[i].fp1 = fp1;
      gz1 = (gzFile) gzopen("chr22_out.fq.gz","wb");
    }    
    else if (output == "fa.gz" || output == "fasta.gz"){
      //FILE *fp1;  error: 'fp1' was not declared in this scope struct_for_threads[i].fp1 = fp1;
      gz1 = (gzFile) gzopen("chr22_out.fa.gz","wb");
    }

    //initialzie values that should be used for each thread
    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;
      struct_for_threads[i].output = output;
      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      struct_for_threads[i].cov = coverage;
      if (output == "fq" || output == "fastq" || output == "fa" || output == "fasta"){
        struct_for_threads[i].fp1 = fp1;
        struct_for_threads[i].gz1 = NULL;}
      else if (output == "fq.gz" || output == "fastq.gz" || output == "fa.gz" || output == "fasta.gz"){
        struct_for_threads[i].fp1 = NULL;
        struct_for_threads[i].gz1 = (gzFile) gz1;}
      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      
      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    //free(catString);
    
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
    if (output == "fq" || output == "fastq" || output == "fa" || output == "fasta"){fclose(fp1);}
    else if (output == "fq.gz" || output == "fastq.gz" || output == "fa.gz" || output == "fasta.gz"){gzclose(gz1);}
    
    
    //for(int i=0;i<nthreads;i++){delete struct_for_threads[i].fqresult_r1 -> s;} //ERROR SUMMARY: 9 errors from 5 contexts (suppressed: 0 from 0 )
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      //ks_release(struct_for_threads[i].fqresult_r1);
      delete struct_for_threads[i].fqresult_r1;
      //free(struct_for_threads[i].fqresult_r1); // ERROR SUMMARY: 8 errors from 4 contexts (suppressed: 0 from 0)
    }
 
    free(genome_data);
  }
  return NULL;
}

// ------------------------------ //
int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // chr1_2.fa  hg19canon.fa  chr1_12   chr22 chr1_15 chr10_15
  //const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa";
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  fprintf(stderr,"\t-> fasta load \n");
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  
  int seed = 1;
  int threads = 5;
  float cov = 1;
  fprintf(stderr,"\t-> Seed used: %d with %d threads\n",seed,threads);

  //char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  const char* Adapt_flag = "false";
  const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
  const char* output = "fa";

  Create_se_threads(seq_ref,threads,seed,cov,output,Adapt_flag,Adapter_1);

  // free the calloc memory from fai_read
  //free(seq_ref);
  fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
  float freq_1 = (float) C_to_T_1.load()/(float) C_total_1.load();
  float freq_2 = (float) C_to_T_2.load()/(float) C_total_2.load();
  float freq_3 = (float) C_to_T_3.load()/(float) C_total_3.load();
  float freq_4 = (float) G_to_A_1.load()/(float) G_total_1.load();
  float freq_5 = (float) G_to_A_2.load()/(float) G_total_2.load();
  float freq_6 = (float) G_to_A_3.load()/(float) G_total_3.load();
  fprintf(stderr,"\t-> Deamination frequency of C-T, pos 1 %.5f , pos 2 %.5f , pos 3 %.5f \n",freq_1,freq_2,freq_3);
  fprintf(stderr,"\t-> Deamination frequency of G-A, pos -1 %.5f , pos -2 %.5f , pos -3 %.5f \n",freq_4,freq_5,freq_6);
}

//SimBriggsModel(seqmod, frag, L, 0.024, 0.36, 0.68, 0.0097);
// g++ SimulAncient_func.cpp atomic.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
//cat test.fq | grep 'T0' | grep 'chr20' | wc -l
//valgrind --tool=memcheck --leak-check=full ./a.out

//awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' chr22_out.fq
// grep 