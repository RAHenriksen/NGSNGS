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

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total){
  // Creates a header for the bamfile. The header is initialized before the function is called //

  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
    
  // Creating header information
  char *name_len_char =(char*) malloc(1024);
  for(int i=0;i<chr_total;i++){
    const char *name = faidx_iseq(seq_ref,i);
    int name_len =  faidx_seq_len(seq_ref,name);
    snprintf(name_len_char,1024,"%d",name_len);
    fprintf(stderr,"ref:%d %d %s\n",i,name,name_len_char);
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
  }
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);
}

std::atomic<float> current_cov_atom(0.0);
std::atomic<int> size_data(0);
std::atomic<int> D_total(0);
std::atomic<int> nread_total(0);

// ---------------------- SINGLE-END ---------------------- //

struct Parsarg_for_Fabam_se_thread{
  // The actual concatenated genome
  char *genome; 
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;

  // File ouput
  samFile *outfile;
  sam_hdr_t *header;

  // Resulting simulated reads ouput
  kstring_t *readid;
  kstring_t *qualstring;
  kstring_t *seq;
  kstring_t *strand;

  const char* Adapter_flag;
  const char* Adapter_1;

  // error profiles
  std::discrete_distribution<> *Ill_err;
  int* sizearray;
  std::discrete_distribution<> *SizeDist;
  const char* read_err_1;

  // simulation parameters
  float current_cov;
  int cov_size;
  int threadseed;
  float cov;
};

void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fabam_se_thread *struct_obj = (Parsarg_for_Fabam_se_thread*) arg;
  //std::cout << "se run "<< std::endl;
  
  // creating random objects for all distributions.
  int seed = struct_obj->threadseed+struct_obj->threadno;
  std::random_device rd;
  std::default_random_engine gen(seed);//gen(struct_obj->threadseed); //struct_obj->seed+struct_obj->threadno

  // -------------------------- // 

  double** R1_2Darray = create2DArray(struct_obj->read_err_1,8,600);
  std::discrete_distribution<> Qualdistr1[600];
  Qual_dist(R1_2Darray,Qualdistr1,600);

  //Free each sub-array
  for(int i = 0; i < 600; ++i) {delete[] R1_2Darray[i];}
  //Free the array of pointers
  delete[] R1_2Darray;

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
    int fraglength = (int) struct_obj->sizearray[struct_obj->SizeDist[1](gen)]; //150; //no larger than 70 due to the error profile which is 280 lines 70 lines for each nt
    
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
      
      int seqlen = strlen(seqmod);

      SimBriggsModel(seqmod, seqmod2, fraglength, 0.024, 0.36, 0.68, 0.0097,std::time(nullptr));

      int strand = rand() % 2;

      if (strand == 0){
        //DNA_complement(seqmod2);
        //reverseChar(seqmod2);
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod2);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          Ill_err(read,struct_obj->Ill_err,gen);
          strncpy(readadapt, read, 150);
          Bam_baseQ(readadapt,qual,Qualdistr1,gen);
          
          ksprintf(struct_obj->readid,"@T%d_RID%d-_S%d_%s:%d-%d_length:%d.",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength);
          ksprintf(struct_obj->qualstring,"%s.",qual);
          ksprintf(struct_obj->seq,"%s.",readadapt);
          ksprintf(struct_obj->strand,"%s.","16");

        }
        else if (struct_obj->Adapter_flag == "false"){
          Ill_err(seqmod2,struct_obj->Ill_err,gen);
          Bam_baseQ(seqmod2,qual,Qualdistr1,gen);

          ksprintf(struct_obj->readid,"@T%d_RID%d-_S%d_%s:%d-%d_length:%d.",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength);
          ksprintf(struct_obj->qualstring,"%s.",qual);
          ksprintf(struct_obj->seq,"%s.",seqmod2);
          ksprintf(struct_obj->strand,"%s.","16");
        }
      }
      else if (strand == 1){
        if (struct_obj->Adapter_flag == "true"){
          strcpy(read, seqmod2);
          strcat(read,struct_obj->Adapter_1);
          //std::cout << "read " << read << std::endl;
          Ill_err(read,struct_obj->Ill_err,gen);
          strncpy(readadapt, read, 150);
          Bam_baseQ(readadapt,qual,Qualdistr1,gen);

          ksprintf(struct_obj->readid,"@T%d_RID%d+_S%d_%s:%d-%d_length:%d.",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength);
          ksprintf(struct_obj->qualstring,"%s.",qual);
          ksprintf(struct_obj->seq,"%s.",readadapt);
          ksprintf(struct_obj->strand,"%s.","0");
        }
        else if (struct_obj->Adapter_flag == "false"){
          Ill_err(seqmod2,struct_obj->Ill_err,gen);
          Bam_baseQ(seqmod2,qual,Qualdistr1,gen);

          ksprintf(struct_obj->readid,"@T%d_RID%d+_S%d_%s:%d-%d_length:%d.",struct_obj->threadno, rand_id,seed,
          struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
          fraglength);
          ksprintf(struct_obj->qualstring,"%s.",qual);
          ksprintf(struct_obj->seq,"%s.",seqmod2);
          ksprintf(struct_obj->strand,"%s.","0");
        }
      }        

      
      bam1_t *bam_file_chr = bam_init1();

      char *token_name, *token_seq, *token_qual, *token_strand;
      char *chr_name;
      char *save_name_ptr, *save_seq_ptr, *save_qual_ptr, *save_strand_ptr;
      char *save_qname_ptr;
      char *qname = (char*) malloc(1024); 
      char *tid_1,*tid_2,*tid_3; 

      token_name = strtok_r(struct_obj->readid->s, ".", &save_name_ptr);
      token_seq = strtok_r(struct_obj->seq->s, ".", &save_seq_ptr);
      token_qual = strtok_r(struct_obj->qualstring->s, ".", &save_qual_ptr);
      token_strand = strtok_r(struct_obj->strand->s, ".", &save_strand_ptr);

      //fprintf(stderr,"read id %s \n",token_name);
      //fprintf(stderr,"seq %s \n",token_seq);
      //fprintf(stderr,"qual %s \n",token_qual);
      
      strcpy(qname, token_name);
      hts_pos_t min_beg, max_end, insert;
      size_t l_qname = strlen(qname);
      uint16_t flag = atoi(token_strand);
      //fprintf(stderr,"%d",flag);
      //fprintf(stderr,"1 %s\n",qname);
      int32_t tid = atoi(strtok_r(qname, ":", &save_qname_ptr));
      //fprintf(stderr,"2 %s\n",qname);
      min_beg = atoi(strtok_r(NULL, "-", &save_qname_ptr))-1; //atoi(strtok_r(NULL, "-", &save_qname_ptr)); 
      //fprintf(stderr,"3 %s\n",qname);
      //fprintf(stderr,"4 %d\n",min_beg);
      tid_1 = qname;
      //fprintf(stderr,"5 %s\n",tid_1);
      while ((tid_2 = strtok_r(NULL, "_", &tid_1))){tid_3 = tid_2;}
      //fprintf(stderr,"5 %s\n",tid_3);
      uint8_t mapq = 60;
      size_t l_aux = 0; // auxiliary field for supp data etc?? 
      //break;
      if (struct_obj->Adapter_flag == "true"){
        char* len_ID = (char*) malloc(1024);
        char* save_len_ptr;

        len_ID = strtok_r(NULL, "", &save_qname_ptr);
        fprintf(stderr," sequence length %s \n",len_ID);
        strtok_r(len_ID, ":", &save_len_ptr);
        int seq_len = atoi(strtok_r(NULL, "", &save_len_ptr));
        uint32_t cigar_bitstring, cigar_bit_soft;
        if (seq_len < 150){
          cigar_bitstring = bam_cigar_gen(seq_len, BAM_CMATCH); // basically does op_len<<BAM_CIGAR_SHIFT|BAM_CMATCH;
          cigar_bit_soft = bam_cigar_gen(strlen(token_seq)-seq_len, BAM_CSOFT_CLIP);
          
        }
        else{
          cigar_bitstring = bam_cigar_gen(150, BAM_CMATCH); // basically does op_len<<BAM_CIGAR_SHIFT|BAM_CMATCH;
          cigar_bit_soft = bam_cigar_gen(0, BAM_CSOFT_CLIP);
        }
        uint32_t cigar_arr[] = {cigar_bitstring,cigar_bit_soft};
        
        size_t n_cigar = 2; // Number of cigar operations, 1 since we only have matches
         //converting uint32_t {aka unsigned int} to const uint32_t* 
        const uint32_t *cigar = cigar_arr;

        pthread_mutex_lock(&Fq_write_mutex);
        bam_set1(bam_file_chr,strlen(token_name),token_name,flag,chr_idx,min_beg,mapq,n_cigar,cigar,-1,-1,0,strlen(token_seq),token_seq,token_qual,l_aux);
        sam_write1(struct_obj->outfile,struct_obj->header,bam_file_chr);
        pthread_mutex_unlock(&Fq_write_mutex);
      }
      else if (struct_obj->Adapter_flag == "false"){
        size_t n_cigar = 1; // Number of cigar operations, 1 since we only have matches

        uint32_t cigar_bitstring = bam_cigar_gen(strlen(token_seq), BAM_CMATCH); // basically does op_len<<BAM_CIGAR_SHIFT|BAM_CMATCH;
        uint32_t cigar_arr[] = {cigar_bitstring}; //converting uint32_t {aka unsigned int} to const uint32_t* 
        const uint32_t *cigar = cigar_arr;

        pthread_mutex_lock(&Fq_write_mutex);
        bam_set1(bam_file_chr,strlen(token_name),token_name,flag,chr_idx,min_beg,mapq,n_cigar,cigar,-1,-1,0,strlen(token_seq),token_seq,token_qual,l_aux);
        sam_write1(struct_obj->outfile,struct_obj->header,bam_file_chr);
        pthread_mutex_unlock(&Fq_write_mutex);
      }
      
      //extract next tokes
      token_name = strtok_r(NULL, ".", &save_name_ptr);
      token_seq = strtok_r(NULL, ".", &save_seq_ptr);
      token_qual = strtok_r(NULL, ".", &save_qual_ptr);
      token_strand = strtok_r(NULL, ".", &save_strand_ptr);
        //token_qual = strtok_r(NULL, "\n", &save_qual_ptr);
        //fprintf(stderr,"%s",token_qual);
      struct_obj->readid->l =0;
      struct_obj->seq->l =0;
      struct_obj->qualstring->l =0;
      struct_obj->strand->l =0;

      nread++;
      //fprintf(stderr,"Number of reads %d \n",nread);
    }

    memset(qual, 0, sizeof(qual));
    memset(readadapt, 0, sizeof seqmod);
    memset(seqmod, 0, sizeof seqmod);
    memset(seqmod2, 0, sizeof seqmod2);

    current_cov_atom = (float) D_total / genome_len;

    struct_obj->current_cov = current_cov_atom; //why do we need this
    struct_obj->cov_size = size_data;

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
  
  //std::cout << "thread done" << std::endl;
  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, float coverage,const char* Adapt_flag,const char* Adapter_1){
  
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];

  //Creates a pointer to allocated memomry for the format
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  if (genome_data != NULL){
    fprintf(stderr,"\t-> Full genome function run!\n");
    fprintf(stderr,"\t-> Full genome size %lu \n",strlen(genome_data));
 
    Parsarg_for_Fabam_se_thread struct_for_threads[nthreads];

    const char *outfile_nam = "Test_pe3.bam";
    samFile *outfile = NULL;
    if ((outfile = sam_open_format(outfile_nam, "wb", fmt_hts)) == 0) {
      fprintf(stderr,"Error opening file for writing\n");
      exit(0);
    }

    // creates a pointer to generated header
    sam_hdr_t *header = sam_hdr_init();
    // add info to the header
    Header_func(fmt_hts,outfile_nam,outfile,header,seq_ref,chr_total);

    // --------------- ERROR PROFILES --------------- //
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

    // --------------- INITIALIZE THREADS --------------- //
    //initialzie values that should be used for each thread

    for (int i = 0; i < nthreads; i++){
      
      struct_for_threads[i].readid =new kstring_t;
      struct_for_threads[i].readid -> l = 0;
      struct_for_threads[i].readid -> m = 0;
      struct_for_threads[i].readid -> s = NULL;

      struct_for_threads[i].qualstring =new kstring_t;
      struct_for_threads[i].qualstring -> l = 0;
      struct_for_threads[i].qualstring -> m = 0;
      struct_for_threads[i].qualstring -> s = NULL;

      struct_for_threads[i].seq =new kstring_t;
      struct_for_threads[i].seq -> l = 0;
      struct_for_threads[i].seq -> m = 0;
      struct_for_threads[i].seq -> s = NULL;

      struct_for_threads[i].strand =new kstring_t;
      struct_for_threads[i].strand -> l = 0;
      struct_for_threads[i].strand -> m = 0;
      struct_for_threads[i].strand -> s = NULL;

      struct_for_threads[i].outfile = outfile;
      struct_for_threads[i].header = header;

      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;

      struct_for_threads[i].Ill_err = Error;
      struct_for_threads[i].sizearray = sizearray;
      struct_for_threads[i].SizeDist = SizeDist;
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      
      struct_for_threads[i].cov = coverage;
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
    
    sam_hdr_destroy(header);
    sam_close(outfile);
    
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].readid -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      delete struct_for_threads[i].readid;
      free(struct_for_threads[i].qualstring -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      delete struct_for_threads[i].qualstring;
      free(struct_for_threads[i].seq -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      delete struct_for_threads[i].seq;
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
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  fprintf(stderr,"\t-> fasta load \n");
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
  
  int seed = 1;
  int threads = 5;
  float cov = 0.5;
  fprintf(stderr,"\t-> Seed used: %d with %d threads\n",seed,threads);

  //char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  const char* Adapt_flag = "false";
  const char* Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";

  Create_se_threads(seq_ref,threads,seed,cov,Adapt_flag,Adapter_1);

  // free the calloc memory from fai_read
  //free(seq_ref);
  fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
}

//SimBriggsModel(seqmod, frag, L, 0.024, 0.36, 0.68, 0.0097);save_len_ptr
// g++ SimulAncient_func.cpp atomic_bam.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
//cat chr22_out.fq | grep '@' | cut -d_ -f4 | sort | uniq -d | wc -l
//cat test.fq | grep 'T0' | grep 'chr20' | wc -l
//valgrind --tool=memcheck --leak-check=full ./a.out

//awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' chr22_out.fq
// grep 
