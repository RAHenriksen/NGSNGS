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
  kstring_t *fqresult_r1;
  faidx_t *seq_ref;
  const char* Ill_err;
  const char* read_err_1;
  const char* read_err_2;
};

void faBam(const char* fastafile,const char* Seq_type,double cov=1.0,const char* r1_profile = "Qual_profiles/Freq_R1.txt",
          const char* r2_profile = "Qual_profiles/Freq_R1.txt"){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
  double init = 1.0;
  int chr_no = 0;
  
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  // ---------------
  // creating the error profiles
  // more generalized way identify number of lines for read profile
  std::ifstream file(r1_profile);
  int Line_no = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();
  size_t lib_size = Line_no/4; 

  // loading in error profiles for nt substitions and read qual
  double** Error_2darray = create2DArray("/home/wql443/WP1/Ill_err.txt",4,280);
  double** R1_2Darray = create2DArray(r1_profile,8,Line_no);
  double** R2_2Darray = create2DArray(r2_profile,8,Line_no);

  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);
  std::discrete_distribution<> Qualdistr2[Line_no];
  Qual_dist(R2_2Darray,Qualdistr2,Line_no);
  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  char qual[1024] = "";
  kstring_t kstr1; 
  kstring_t kstr2;
  
  // ---------------
  //create the file for saving

  //Creates a pointer to allocated memomry for the format??
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
    
  // Creating sam_header
  char *name_len_char =(char*) malloc(1024);
  for(int i=0;i<faidx_nseq(seq_ref);i++){
    const char *name = faidx_iseq(seq_ref,i);
    int name_len =  faidx_seq_len(seq_ref,name);
    // skal være c string i array så vi konvertere int om til char
    snprintf(name_len_char,1024,"%d",name_len);
    fprintf(stderr,"ref:%d %d %s\n",i,name,name_len_char);
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
  }
  
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);
  
  // alignment delen, gemt i bam_1 type for at representere hver linje som 1 alignment
  bam1_t *bam_file = bam_init1(); //initialisere bam1_t (type) til hukommelsen!

  // ---------------
  // extracting the sequences
  while (chr_no < faidx_nseq(seq_ref)){  
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);

    int start_pos = 1;
    int end_pos = name_len; //30001000
    char seqmod[1024];

    while(start_pos <= end_pos){

      int readlength = drand48()*(70.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;

      //extracts the sequence
      strncpy(seqmod,data+start_pos,readlength);

      char Qname[96]; //read_id
      snprintf(Qname,96,"%s:%d-%d",name,start_pos,stop);

      char * pch;
      pch = strchr(seqmod,'N');

      if (pch != NULL){start_pos += readlength + 1;}
      else {
        Ill_err(seqmod,Error,gen);
        Bam_baseQ(seqmod,qual,Qualdistr1,gen);

        //ensures proper bam format with 11 mandatory fields mandatory fields
        // QNAME = char Qname[96] 
        int Flag = 4; // 4 for unmapped
        int RNAME = chr_no; // Reference sequence name, chr_no takes chr from sam_hdr_t
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
        bam_set1(bam_file,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,chr_no,Pnext-1,Tlen,readlength,seqmod,qual,0);
        // -1 for the different positions is because they are made 1 - based with the bam_set
        
        sam_write1(outfile,header,bam_file);
        memset(qual, 0, sizeof(qual));
      }
      start_pos += readlength + 1;
    }
    chr_no++;
  }
  sam_hdr_destroy(header);
  sam_close(outfile);
  return; 
}


int main(int argc,char **argv){
  clock_t tStart = clock();
  // /willerslev/users-shared/science-snm-willerslev-wql443/reference_files/Human/hg19canon.fa
  // "/home/wql443/scratch/reference_genome/hg19/chr22.fa"
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr2122.fa";
  int seed = -1;
  if (std::strcmp(argv[1],"bam")==0){
    if(argc>3){seed = atoi(argv[3]);}
    else{seed = time(NULL);} 
    fprintf(stderr,"seed is: %d\n",seed);
    srand48(seed);
    faBam(fastafile,"pe");
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}

// g++ SimulAncient_func.cpp FaBam_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -Wall
