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
  kstring_t *name;
  kstring_t *qualstring;
  kstring_t *seq;
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
      int len = start_pos+readlength;
      ksprintf(struct_obj->name,"%s.",Qname);
      ksprintf(struct_obj->qualstring,"%s.",qual);
      ksprintf(struct_obj->seq,"%s.",seqmod);
      //ksprintf(struct_obj->seq,"%s,%s,%s.",Qname,seqmod,qual);
      //ksprintf(struct_obj->seq,"%s_",Qname,"%s.",seqmod);
      
    }
    memset(seqmod, 0, sizeof(seqmod));
    memset(qual, 0, sizeof(qual));
    start_pos += readlength + 1;
  }
}

/*
int main(int argc,char **argv){
  //Loading in an creating my objects for the sequence files.
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr20_22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = faidx_nseq(seq_ref);
  
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
  }
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);

  // Initializing the alignment information in a bam_1 type to memory, with each line representing one alignment 
  //bam1_t *bam_file = bam_init1();

  //initialize mutex
  pthread_mutex_init(&data_mutex,NULL);
  int nthreads=chr_no;
  pthread_t mythreads[nthreads];
  Parsarg_for_Fabam_thread struct_for_threads[nthreads];
  
  //initialzie values that should be used for each thread
  
  //bam1_t *bam_file_chr = ;  // Saves 25 lines each with same final output from last chromosome it iterates throug
  for(int i=0;i<nthreads;i++){
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
    struct_for_threads[i].read_err_2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";
  }
  //struct_for_threads[i].bam_file = bam_file_chr(); 
  //    struct_for_threads[i].out_bam_name = outfile;  struct_for_threads[i].bam_head = header;
  //launch all worker threads
  for(int i=0;i<nthreads;i++){
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_create(&mythreads[i],&attr,FaBam_thread_run,&struct_for_threads[i]);
  }

  //right now 22 thread are running at the same time
  //now join, this means waiting for all worker threads to finish
  for(int i=0;i<nthreads;i++){
    pthread_join(mythreads[i],NULL);
  }

  pthread_mutex_destroy(&data_mutex);  
  //now all are done and we only write in main program.
  for(int i=0;i<nthreads;i++){
    bam1_t *bam_file_chr = bam_init1();
    //str
    char * pch_seq;
    char * pch_name;
    char * pch_qual;
    pch_name = strtok(struct_for_threads[i].name->s,".");
    pch_seq = strtok(struct_for_threads[i].seq->s,".");
    pch_qual = strtok(struct_for_threads[i].qualstring->s,".");

    std::cout << "-------" << std::endl;
    std::cout << pch_seq << std::endl;
    //pch_seq = strtok(struct_for_threads[i].seq->,".");
    //pch_qual = strtok(struct_for_threads[i].qual->s,".");
    bam_set1(bam_file_chr,strlen(pch_name),pch_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(pch_seq),pch_seq,pch_qual,0);
    sam_write1(outfile,header,bam_file_chr);
    / *
    while (pch_qual != NULL & pch_name!= NULL)
    {
      //printf("%s\n",pch_name);
      //bam_set1(bam_file_chr,strlen(pch_name),pch_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(pch_seq),pch_seq,pch_seq,0);
      bam_set1(bam_file_chr,strlen(pch_name),pch_name,4,-1,-1,0,0,NULL,-1,-1,0,strlen(pch_qual),pch_qual,pch_qual,0);
      sam_write1(outfile,header,bam_file_chr);
      pch_qual = strtok (NULL, ".");
      //pch_seq = strtok (NULL, ".");
      pch_name = strtok (NULL, ".");

    }
    * /
  }
  sam_hdr_destroy(header);
  sam_close(outfile);
} 
*/
// g++ SimulAncient_func.cpp FaBam_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall

/*
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
      struct_obj->x->name=Qname;
      struct_obj->x->qual=qual;
      struct_obj->x->seq = seqmod;

          sam_write1(outfile,header,bam_file);
  sam_hdr_destroy(header);
  sam_close(outfile);
  return; */

/*
int main ()
{
  char str[] ="AAAAGAGAGAGAGAGA.>>>>>>>>>>>>>>";
  char str2[] ="CCCCCCCCCCCCCCC.XXXXXXXXXXXXX.AAAAAAA";
  char * pch;
  pch = strtok (str," .");
  char * pch2;
  pch2 = strtok (str2," .");
  char *str3 = "CCCCCCCCCCCCCCC.XXXXXXXXXXXXX";
  char *last = strrchr(str3, '.');
  printf("Last token: '%s'\n", last);
  printf("Last token: '%s'\n", last+1);
  printf("Last token: '%s'\n", last+2);
  while (pch != NULL)
  {
    //printf ("%s\n",pch);
    std::cout << pch << std::endl;
    pch = strtok (NULL, " .");
    //std::cout << pch2 << std::endl;
    //pch2 = strtok (NULL, " .");
  }
  return 0;
}*/

/*
int main()
{
    kstring_t kstr1;
    kstr1.s = NULL; kstr1.l = 0; kstr1.m = 0;
    char input[1024] = "READID1_AGAGAG_BBBBBBB.READID2_GGGGGGGG_FFFFFFFF.READID3_CCC_III";
    ksprintf(&kstr1,"%s",input);
    char *tokens1;
    tokens1 = strtok(kstr1.s,".");
    while (tokens1 != NULL)
    {
      std::cout << "while 1 " << tokens1 << std::endl;
      tokens1 = strtok(NULL, ".");
    }
}*/

int main()
{
    kstring_t kstr1;kstring_t kstr2;kstring_t kstr3;
    kstr1.s = NULL; kstr1.l = 0; kstr1.m = 0;
    kstr2.s = NULL; kstr2.l = 0; kstr2.m = 0;
    kstr3.s = NULL; kstr3.l = 0; kstr3.m = 0;
    char name[1024] = "ReadID1.ReadID2.ReadID3";
    char seq[1024] = "AAAA.GGGG.CCCC";
    char qual[1024] = "BBBB.FFFF.IIII";

    ksprintf(&kstr1,"%s",name);
    ksprintf(&kstr2,"%s",seq);
    ksprintf(&kstr3,"%s",qual);

    char *token1; char *token2;char *token3;
    char *save_ptr1, *save_ptr2, *save_ptr3;

    token1 = strtok_r(name, ".", &save_ptr1);
    token2 = strtok_r(seq, ".", &save_ptr2);
    token3 = strtok_r(qual, ".", &save_ptr3);

    while(token1 != NULL && token2 != NULL && token3 != NULL) {
      std::cout << token1 << std::endl;
      std::cout << token2 << std::endl;
      std::cout << token3 << std::endl;
      // get next tokens
      token1 = strtok_r(NULL, ".", &save_ptr1);
      token2 = strtok_r(NULL, ".", &save_ptr2);
      token3 = strtok_r(NULL, ".", &save_ptr3);
      std::cout << "while loop end "<< std::endl;
    }
}
/*
    char name[1024] = "READID1.READID2.READID3";
    char seq[1024] = "AAAA_GGGGGG_TTTTTTTT";
    char qual[1024] = "BBBB,FFFFFF,>>>>>>>>";
    ksprintf(&kstr1,"%s",name);
    ksprintf(&kstr2,"%s",seq);
    ksprintf(&kstr3,"%s",qual);
    std::cout << kstr1.s << std::endl;
    char *tokens1; char *tokens2; char *tokens3;
    tokens1 = strtok(kstr1.s,".");
    std::cout << tokens1 << std::endl;
    tokens2 = strtok(kstr2.s,"_");
    tokens3 = strtok(kstr3.s,",");
    char * test;
    while (tokens1 != NULL)
    {
      std::cout << "while 1 " << tokens1 << std::endl;
      tokens1 = strtok(NULL, ".");
    }
*/