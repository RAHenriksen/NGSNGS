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

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>          // std::mutex mtx;

#include "SimulAncient_func.h"

pthread_mutex_t data_mutex;

struct Parsarg_for_Fabam_thread{
  int chr_idx;
  faidx_t *seq_ref;
  bam1_t *bam_format;
  samFile *bam_out;
  sam_hdr_t *bam_head;
  const char* Ill_err;
  const char* read_err_1;
  const char* read_err_2;
};
//  sam_hdr_t *bam_head;  samFile *out_bam_name;
void* FaBam_thread_run(void *arg){
  
  Parsarg_for_Fabam_thread *struct_obj = (Parsarg_for_Fabam_thread*) arg;

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
  int end_pos = chr_len;

  char seqmod[1024] = {0};
  char qual[1024] = "";

  while(start_pos <= end_pos){
    int readlength = drand48()*(70.0-30.0)+30.0;
    int stop = start_pos+(int) readlength;

    //extracts the sequence
    strncpy(seqmod,data+start_pos,readlength);

    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += readlength + 1;}
    else{ 
      char Qname[96]; //read_id
      snprintf(Qname,96,"%s:%d-%d",chr_name,start_pos,stop);

      Ill_err(seqmod,Error,gen);
      Bam_baseQ(seqmod,qual,Qualdistr1,gen);

      //ensures proper bam format with 11 mandatory fields mandatory fields
      // QNAME = char Qname[96] 
      int Flag = 4; // 4 for unmapped
      int RNAME = idx; // Reference sequence name, chr_no takes chr from sam_hdr_t
      // POS = int start_pos -1 as it is already 1-based from faidx but bam_set1 converts it further to 1-based
      int mapQ = 255; // 255 if unavailable du to no mapping
      const uint32_t *cigar = NULL; // cigar string, NULL if unavailable - But how do we then add actual cigar info?
      // RNEXT = mtid -> chr_no
      int Pnext = 0; // position for next mate -> mpos
      int Tlen = 0; // template length -> isize
      // SEQ = sequence
      // QUAL = quality
      int no_cigar = 0; //number of cigar operations 
      
      //bam, lqname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux
      //struct_obj->bam_line

      pthread_mutex_lock(&data_mutex);
      bam_set1(struct_obj->bam_format,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,idx,Pnext-1,Tlen,readlength,seqmod,qual,0);
      sam_write1(struct_obj->bam_out,struct_obj->bam_head,struct_obj->bam_format);
      pthread_mutex_unlock(&data_mutex);
        
      // -1 for the different positions is because they are made 1 - based with the bam_set
    }
    memset(qual, 0, sizeof(qual));
    start_pos += readlength + 1;
  }
  pthread_exit(0);
}


int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  // /home/wql443/scratch/reference_genome/hg19/chr2122.fa
  // /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa
  // chr1_2
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr21_22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  std::cout << "chromosome number " << chr_no <<std::endl;
  // Initializing the bam file header

  //Creates a pointer to allocated memomry for the format
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  //wb -> bam , wc -> cram
  char out_mode[5]="wb";
    
  const char *outfile_nam = "Test_pe.bam";
  samFile *outfile = NULL;

  if ((outfile = sam_open_format(outfile_nam, out_mode, fmt_hts)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(0);
  }
  
  // creates a pointer to generated header
  sam_hdr_t *header = sam_hdr_init();
  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
  char *refName=NULL;
    
  // Creating header information
  char *name_len_char =(char*) malloc(1024);
  for(int i=0;i<chr_no;i++){
    const char *name = faidx_iseq(seq_ref,i);
    int name_len =  faidx_seq_len(seq_ref,name);
    snprintf(name_len_char,1024,"%d",name_len);
    fprintf(stderr,"ref:%d %d %s\n",i,name,name_len_char);
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
    std::cout << "done" <<std::endl;
  }
  free(name_len_char);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);

  // Initializing the alignment information in a bam_1 type to memory, with each line representing one alignment 
  //bam1_t *bam_file = bam_init1();

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  int nthreads=chr_no;
  pthread_t mythreads[nthreads];

  Parsarg_for_Fabam_thread struct_for_threads[nthreads]; //should go to the stack

  // creating a pointer to allocated memory by type casting my struct and malloc with the size of all the elements in my struct
  //struct Parsarg_for_Fabam_thread *struct_for_threads = (struct Parsarg_for_Fabam_thread*) malloc (chr_no * sizeof (struct Parsarg_for_Fabam_thread));

  // New is a similar to malloc
  //Parsarg_for_Fabam_thread* struct_for_threads = new Parsarg_for_Fabam_thread[nthreads];

  //initialzie values that should be used for each thread
  bam1_t *bam_file_chr =bam_init1() ;  // Saves 25 lines each with same final output from last chromosome it iterates throug
  
  for(int i=0;i<nthreads;i++){
    struct_for_threads[i].chr_idx = i;
    struct_for_threads[i].seq_ref = seq_ref;
    struct_for_threads[i].bam_format = bam_file_chr; //bam_init() saves final output from all chromosomes -> still 25 lines
    struct_for_threads[i].bam_head=header;
    struct_for_threads[i].bam_out=outfile;
    struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
    struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
    struct_for_threads[i].read_err_2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_create(&mythreads[i],&attr,FaBam_thread_run,&struct_for_threads[i]);
  }
  // free(struct_for_threads);
  //struct_for_threads[i].bam_file = bam_file_chr(); 
  //struct_for_threads[i].out_bam_name = outfile;  struct_for_threads[i].bam_head = header;

  //launch all worker threads

  pthread_mutex_destroy(&data_mutex);  
  //right now 22 thread are running at the same time
  //now join, this means waiting for all worker threads to finish
  for(int i=0;i<nthreads;i++){
    pthread_join(mythreads[i],NULL);
    bam_destroy1(struct_for_threads[i].bam_format);
  }
  //now all are done and we only write in main program.
  //for(int i=0;i<nthreads;i++){
  //  sam_write1(outfile,header,struct_for_threads[i].bam_file);
  //}
  sam_hdr_destroy(header);
  sam_close(outfile);
} 
// g++ SimulAncient_func.cpp FaBam_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall

/*

          sam_write1(outfile,header,bam_file);
  sam_hdr_destroy(header);
  sam_close(outfile);
  return; */