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

// -------------- SINGLE END DATA ------------- //

struct Parsarg_for_Fabam_se_thread{
  int chr_idx;
  faidx_t *seq_ref;
  kstring_t *name;
  kstring_t *qualstring;
  kstring_t *seq;
  const char* Ill_err;
  const char* read_err_1;
};
//  sam_hdr_t *bam_head;  samFile *out_bam_name;
void* FaBam_thread_se_run(void *arg){
  
  Parsarg_for_Fabam_se_thread *struct_obj = (Parsarg_for_Fabam_se_thread*) arg;

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

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  int start_pos = 1;
  int end_pos = chr_len;

  char seqmod[1024] = {0};
  char qual[1024] = "";

  while(start_pos <= end_pos){
    int fraglength = (int) sizearray[SizeDist[1](gen)];
    //int readlength = drand48()*(70.0-30.0)+30.0;
    int stop = start_pos+(int) fraglength;

    // extract the sequences and limits to a maximum size of 150
    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod,data+start_pos,150);
    }
    // case 2
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod,data+start_pos,150);
    }
    // case 3
    else if (fraglength <= 150)
    {
      strncpy(seqmod,data+start_pos,fraglength);
    }
    //extracts the sequence
    //strncpy(seqmod,data+start_pos,readlength);

    //removes NNN
    char * pch;
    pch = strchr(seqmod,'N');
    if (pch != NULL){start_pos += fraglength + 1;}
    else{ 
      char Qname[96]; //read_id
      snprintf(Qname,96,"%s:%d-%d_length:%d",chr_name,start_pos+1,stop,fraglength);
      Ill_err(seqmod,Error,gen);
      Bam_baseQ(seqmod,qual,Qualdistr1,gen);
      /*std::cout << Qname << std::endl;
      std::cout << seqmod << std::endl;
      std::cout << qual << std::endl;*/
      
      //ensures proper bam format with 11 mandatory fields mandatory fields
      // QNAME = char Qname[96] 
 
      
      //bam, lqname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux
      //struct_obj->bam_line
      //bam_set1(struct_obj->bam_format,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,idx,Pnext-1,Tlen,readlength,seqmod,qual,0);
      //sam_write1(struct_obj->bam_out,struct_obj->bam_head,struct_obj->bam_format);
      // -1 for the different positions is because they are made 1 - based with the bam_set
      int len = start_pos+fraglength;
      ksprintf(struct_obj->name,"%s.",Qname);
      ksprintf(struct_obj->qualstring,"%s.",qual);
      ksprintf(struct_obj->seq,"%s.",seqmod);
      //ksprintf(struct_obj->seq,"%s,%s,%s.",Qname,seqmod,qual);
      //ksprintf(struct_obj->seq,"%s_",Qname,"%s.",seqmod);
      
    }
    memset(seqmod, 0, sizeof(seqmod));
    memset(qual, 0, sizeof(qual));
    start_pos += fraglength;
  }
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no,int chr_total){
  int chr_idx = 0;

  // Initializing the bam file header

  //Creates a pointer to allocated memomry for the format
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  
  //wb -> bam , wc -> cram
  char out_mode[5]="wb";
    
  const char *outfile_nam = "Test_pe3.bam";
  samFile *outfile = NULL;

  if ((outfile = sam_open_format(outfile_nam, out_mode, fmt_hts)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(0);
  }
  
  // creates a pointer to generated header
  sam_hdr_t *header = sam_hdr_init();
  // add info to the header
  Header_func(fmt_hts,outfile_nam,outfile,header,seq_ref,chr_total);

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);

  int nthreads = chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fabam_se_thread struct_for_threads[nthreads];
  
  int i = 0; int j = 0; int k = 0;
  while(chr_idx < nthreads){
    //initialize values that should be used for each thread
    for(i; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].name =new kstring_t;
      struct_for_threads[i].name -> l = 0;
      struct_for_threads[i].name -> m = 0;
      struct_for_threads[i].name -> s = NULL;
      struct_for_threads[i].qualstring=new kstring_t;
      struct_for_threads[i].qualstring -> l = 0;
      struct_for_threads[i].qualstring -> m = 0;
      struct_for_threads[i].qualstring -> s = NULL;
      struct_for_threads[i].seq=new kstring_t;
      struct_for_threads[i].seq -> l = 0;
      struct_for_threads[i].seq -> m = 0;
      struct_for_threads[i].seq -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      std::cout << "chr_no" << i << std::endl;
    }
    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      std::cout << "launc for loop " << std::endl;
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,FaBam_thread_se_run,&struct_for_threads[j]);
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
  for(int i=0;i<nthreads;i++){
    bam1_t *bam_file_chr = bam_init1();
    
    char *token_name; char *token_seq;char *token_qual;
    char *save_name_ptr, *save_seq_ptr, *save_qual_ptr;
    token_name = strtok_r(struct_for_threads[i].name->s, ".", &save_name_ptr);
    token_seq = strtok_r(struct_for_threads[i].seq->s, ".", &save_seq_ptr);
    token_qual = strtok_r(struct_for_threads[i].qualstring->s, ".", &save_qual_ptr);
    //std::cout << struct_for_threads[i].name->l << std::endl;
    //std::cout << struct_for_threads[i].seq->l << std::endl;
    
    while(token_name != NULL && token_seq != NULL && token_qual != NULL) {
      //std::cout << token_seq<< std::endl;
      bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(token_seq),token_seq,token_qual,0);
      sam_write1(outfile,header,bam_file_chr);
      //extract next tokes
      token_name = strtok_r(NULL, ".", &save_name_ptr);
      token_seq = strtok_r(NULL, ".", &save_seq_ptr);
      token_qual = strtok_r(NULL, ".", &save_qual_ptr);
      //std::cout << "while loop end "<< std::endl;
    }
  }

  sam_hdr_destroy(header);
  sam_close(outfile);
}

// -------------- PAIRED END DATA ------------- //


struct Parsarg_for_Fabam_pe_thread{
  int chr_idx;
  faidx_t *seq_ref;
  kstring_t *name;
  kstring_t *qualstring_r1;
  kstring_t *qualstring_r2;
  kstring_t *seq_r1;
  kstring_t *seq_r2;
  const char* Ill_err;
  const char* read_err_1;
  const char* read_err_2;
};
//  sam_hdr_t *bam_head;  samFile *out_bam_name;
void* FaBam_thread_pe_run(void *arg){
  
  Parsarg_for_Fabam_pe_thread *struct_obj = (Parsarg_for_Fabam_pe_thread*) arg;

  int idx = struct_obj->chr_idx;
  const char *chr_name = faidx_iseq(struct_obj->seq_ref,idx);
  int chr_len = faidx_seq_len(struct_obj->seq_ref,chr_name);
  
  pthread_mutex_lock(&data_mutex);
  char *data = fai_fetch(struct_obj->seq_ref,chr_name,&chr_len);
  pthread_mutex_unlock(&data_mutex);

  //just load in first error file to count number of lines
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

  // quality error distribution for read 1
  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);

  // quality error distribution for read 2
  std::discrete_distribution<> Qualdistr2[Line_no];
  Qual_dist(R2_2Darray,Qualdistr2,Line_no);
  
  // sequencing error distribution
  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  // Creates the random lengths array and distributions //
  std::ifstream infile("Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();

  std::discrete_distribution<> SizeDist[2]; 
  
  std::ifstream infile2("Size_freq.txt");
  Size_freq_dist(infile2,SizeDist); //creates the distribution of all the frequencies
  infile2.close();
  // ---------------------- //

  int start_pos = 1;
  int end_pos = chr_len;

  char seqmod_r1[1024] = {0};
  char seqmod_r2[1024] = {0};
  char qual_r1[1024] = "";
  char qual_r2[1024] = "";

  while(start_pos <= end_pos){
    int fraglength = (int) sizearray[SizeDist[1](gen)];
    //int readlength = drand48()*(70.0-30.0)+30.0;
    int stop = start_pos+(int) fraglength;

    // extract the sequences and limits to a maximum size of 150
    // case 1
    if (fraglength > 2*150){
      //std::cout << "lolrt" << std::endl;
      strncpy(seqmod_r1,data+start_pos,150);
      strncpy(seqmod_r2,data+stop-150,150);
    }
    // case 2
    else if (150 < fraglength && fraglength < 2*150) //case 2
    {
      strncpy(seqmod_r1,data+start_pos,150);
      strncpy(seqmod_r2,data+stop-150,150);
    }
    // case 3
    else if (fraglength <= 150)
    {
      strncpy(seqmod_r1,data+start_pos,fraglength);
      strncpy(seqmod_r2,data+start_pos,fraglength);
    }
    //extracts the sequence
    //strncpy(seqmod_r1,data+start_pos,readlength);

    //removes NNN
    char * pch;
    pch = strchr(seqmod_r1,'N');
    if (pch != NULL){start_pos += fraglength + 1;}
    else{ 
      char Qname[96]; //read_id
      snprintf(Qname,96,"%s:%d-%d_length:%d",chr_name,start_pos+1,stop,fraglength);
      Ill_err(seqmod_r1,Error,gen);
      Bam_baseQ(seqmod_r1,qual_r1,Qualdistr1,gen);
      
      Ill_err(seqmod_r2,Error,gen);
      Bam_baseQ(seqmod_r2,qual_r2,Qualdistr2,gen);

      ksprintf(struct_obj->name,"%s.",Qname);
      ksprintf(struct_obj->qualstring_r1,"%s.",qual_r1);
      ksprintf(struct_obj->seq_r1,"%s.",seqmod_r1);
      ksprintf(struct_obj->qualstring_r2,"%s.",qual_r2);
      ksprintf(struct_obj->seq_r2,"%s.",seqmod_r2);
    }
    memset(seqmod_r1, 0, sizeof(seqmod_r1));
    memset(seqmod_r2, 0, sizeof(seqmod_r2));
    memset(qual_r1, 0, sizeof(qual_r1));
    memset(qual_r2, 0, sizeof(qual_r2));    
    start_pos += fraglength;
  }
}

void* Create_pe_threads(faidx_t *seq_ref,int thread_no,int chr_total){
  int chr_idx = 0;
  std::cout << "Create_Pe_Trehads" << std::endl;
  // Initializing the bam file header

  //Creates a pointer to allocated memomry for the format
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  
  //wb -> bam , wc -> cram
  char out_mode[5]="wb";
    
  const char *outfile_nam = "Test_PE.bam";
  samFile *outfile = NULL;

  if ((outfile = sam_open_format(outfile_nam, out_mode, fmt_hts)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(0);
  }
  
  // creates a pointer to generated header
  sam_hdr_t *header = sam_hdr_init();
  // add info to the header
  Header_func(fmt_hts,outfile_nam,outfile,header,seq_ref,chr_total);

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);

  int nthreads = chr_total;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fabam_pe_thread struct_for_threads[nthreads];
  
  int i = 0; int j = 0; int k = 0;
  while(chr_idx < nthreads){
    //initialize values that should be used for each thread
    for(i; i < std::min(chr_idx+thread_no,nthreads);i++){
      struct_for_threads[i].chr_idx = i;
      struct_for_threads[i].seq_ref = seq_ref;
      struct_for_threads[i].name =new kstring_t;
      struct_for_threads[i].name -> l = 0;
      struct_for_threads[i].name -> m = 0;
      struct_for_threads[i].name -> s = NULL;
      struct_for_threads[i].qualstring_r1 = new kstring_t;
      struct_for_threads[i].qualstring_r1 -> l = 0;
      struct_for_threads[i].qualstring_r1 -> m = 0;
      struct_for_threads[i].qualstring_r1 -> s = NULL;
      struct_for_threads[i].seq_r1 = new kstring_t;
      struct_for_threads[i].seq_r1 -> l = 0;
      struct_for_threads[i].seq_r1 -> m = 0;
      struct_for_threads[i].seq_r1 -> s = NULL;
      struct_for_threads[i].qualstring_r2 = new kstring_t;
      struct_for_threads[i].qualstring_r2 -> l = 0;
      struct_for_threads[i].qualstring_r2 -> m = 0;
      struct_for_threads[i].qualstring_r2 -> s = NULL;
      struct_for_threads[i].seq_r2 = new kstring_t;
      struct_for_threads[i].seq_r2 -> l = 0;
      struct_for_threads[i].seq_r2 -> m = 0;
      struct_for_threads[i].seq_r2 -> s = NULL;
      struct_for_threads[i].Ill_err = "/home/wql443/WP1/SimulAncient/Qual_profiles/Ill_err.txt";
      struct_for_threads[i].read_err_1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
      struct_for_threads[i].read_err_2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";
      std::cout << "chr_no" << i << std::endl;
    }
    //launch worker threads
    for (j; j < std::min(chr_idx+thread_no,nthreads);j++){
      std::cout << "launc for loop " << std::endl;
      pthread_attr_t attr;
      pthread_attr_init(&attr);
      pthread_create(&mythreads[j],&attr,FaBam_thread_pe_run,&struct_for_threads[j]);
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
  for(int i=0;i<nthreads;i++){
    bam1_t *bam_file_chr = bam_init1();
    
    char *token_name; char *token_seq_r1;char *token_qual_r1;
    char *save_name_ptr; char *save_seq_r1_ptr; char *save_qual_r1_ptr;
    
    char *token_seq_r2; char *token_qual_r2;
    char *save_seq_r2_ptr; char *save_qual_r2_ptr;

    //char *qname;

    token_name = strtok_r(struct_for_threads[i].name->s, ".", &save_name_ptr);
    //qname = strtok_r(struct_for_threads[i].name->s, ".", &save_qname_ptr);
    token_seq_r1 = strtok_r(struct_for_threads[i].seq_r1->s, ".", &save_seq_r1_ptr);
    token_qual_r1 = strtok_r(struct_for_threads[i].qualstring_r1->s, ".", &save_qual_r1_ptr);
    token_seq_r2 = strtok_r(struct_for_threads[i].seq_r2->s, ".", &save_seq_r2_ptr);
    token_qual_r2 = strtok_r(struct_for_threads[i].qualstring_r2->s, ".", &save_qual_r2_ptr);
    
    //std::cout << struct_for_threads[i].name->l << std::endl;
    //std::cout << struct_for_threads[i].seq->l << std::endl;
    char *save_qname_ptr;
    char* qname = (char*) malloc(1024); 
    strcpy(qname, token_name);

    std::cout << "CIFGAR" << std::endl;
    while(token_name != NULL && token_seq_r1 != NULL && token_qual_r1 != NULL) {
      
      //std::cout << token_name<< std::endl;
      //std::cout << strtok(token_name,":") << std::endl;
      //std::cout << strtok(NULL,"-") << std::endl;
      /*int a = strlen(token_seq_r1);
      const uint32_t arr[] = {a};
      const uint32_t *cigar = arr;
      /std::cout << "-----------" << std::endl;
      std::cout << cigar << std::endl;
      std::cout << &cigar << std::endl;
      std::cout << *cigar << std::endl;
      std::cout << a << " " << token_seq_r1 << std::endl;*/

      hts_pos_t min_beg, max_end;
      
      size_t l_qname = strlen(qname);
      uint16_t flag = 4;
      //std::cout << "qname before " << qname << std::endl;
      int32_t tid = atoi(strtok_r(qname, ":", &save_qname_ptr));
      /*std::cout << "qname after " << qname << std::endl;
      std::cout << strtok_r(NULL, "-", &save_qname_ptr) << std::endl;
      std::cout << "qname after 2 " << qname << std::endl;*/
      min_beg = 100; //atoi(strtok_r(NULL, "-", &save_qname_ptr)); 
      uint8_t mapq = 60;
      size_t n_cigar = strlen(token_seq_r1);
      uint32_t arr[] = {0};
      //memset(arr, 0, sizeof(arr)*n_cigar);
      const uint32_t *cigar = arr;
      size_t l_aux = 0; // auxiliary field for supp data etc?? 
      bam_set1(bam_file_chr,strlen(token_name),token_name,flag,-1,min_beg,mapq,n_cigar,cigar,tid,-1,0,strlen(token_seq_r1),token_seq_r1,token_qual_r1,l_aux);
      sam_write1(outfile,header,bam_file_chr);
      
      /*bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(token_seq_r1),token_seq_r1,token_qual_r1,0);
      sam_write1(outfile,header,bam_file_chr);
      bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(token_seq_r1),token_seq_r1,token_qual_r1,0);
      sam_write1(outfile,header,bam_file_chr);*/
      
      // MEN DISSE TOKENS OVERSKRIVER DE GAMLE!
      //extract next tokes
      token_name = strtok_r(NULL, ".", &save_name_ptr);
      token_seq_r1 = strtok_r(NULL, ".", &save_seq_r1_ptr);
      token_qual_r1 = strtok_r(NULL, ".", &save_qual_r1_ptr);
      token_seq_r2 = strtok_r(NULL, ".", &save_seq_r2_ptr);
      token_qual_r2 = strtok_r(NULL, ".", &save_qual_r2_ptr);
      //std::cout << "while loop end "<< std::endl;
    }
  }

  sam_hdr_destroy(header);
  sam_close(outfile);
}


// ------------- MAIN ------------- // 
int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  //chr12_15
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_total = faidx_nseq(seq_ref);
  
  // Creates the threads and writes to the files
  //Create_se_threads(seq_ref,chr_total,chr_total);
  Create_pe_threads(seq_ref,chr_total,chr_total);

} 

// g++ SimulAncient_func.cpp FaBam_kstrings.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl


//bam, lqname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux
//bam_set1(struct_obj->bam_format,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,idx,Pnext-1,Tlen,readlength,seqmod,qual,0);

/*
while(token_name != NULL && token_seq_r1 != NULL && token_qual_r1 != NULL) {
      std::cout << "while loop " << std::endl;
      int a = strlen(token_seq_r1);
      const uint32_t arr[] = {a};
      const uint32_t *cigar = arr;
      //std::cout << token_seq<< std::endl;
      //bam, lqname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux
      //bam_set1(struct_obj->bam_format,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,idx,Pnext-1,Tlen,readlength,seqmod,qual,0);
      //bam_set1(bam_file_chr,strlen(token_name),token_name,99,-1,-1,60,strlen(token_seq_r1),&cigar,-1,-1,0,strlen(token_seq_r1),token_seq_r1,token_qual_r1,0);
      //bam_set1(bam_file_chr,strlen(token_name),token_name,147,-1,-1,60,strlen(token_seq_r2),cigar,-1,-1,0,strlen(token_seq_r2),token_seq_r2,token_qual_r2,0);
      bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,1,NULL,-1,-1,0,strlen(token_seq_r1),token_seq_r1,token_qual_r1,0);
      //bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,1,NULL,-1,-1,0,strlen(token_seq_r2),token_seq_r2,token_qual_r2,0);
      //bam_set1(bam_file_chr,strlen(token_name),token_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(token_seq),token_seq,token_qual,0);
      sam_write1(outfile,header,bam_file_chr);
      
      // if the flags are 99 and 147 then i need the CIGAR STRING

      //extract next tokes
      token_name = strtok_r(NULL, ".", &save_name_ptr);
      token_seq_r1 = strtok_r(NULL, ".", &save_seq_r1_ptr);
      token_qual_r1 = strtok_r(NULL, ".", &save_qual_r1_ptr);
      token_seq_r2 = strtok_r(NULL, ".", &save_seq_r2_ptr);
      token_qual_r2 = strtok_r(NULL, ".", &save_qual_r2_ptr);
    }*/