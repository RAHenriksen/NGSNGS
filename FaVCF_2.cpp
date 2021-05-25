#include <cstdio>
#include <cassert>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <typeinfo>
#include <random>
#include <iterator>
#include <cmath>
#include <chrono>
#include <time.h>
#include <algorithm>

void faVCF(const char* fastafile,FILE *fp){
  //creates fasta structure
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;
  
  // VCF STRUCTURE
  htsFile *test_bcf = NULL;
  bcf_hdr_t *test_header = NULL;

  test_bcf = bcf_open("/home/wql443/WP1/data/Chr14_subset.vcf", "r");
  test_header = bcf_hdr_read(test_bcf);
  bcf1_t *test_record = bcf_init();
  int vcfstruct = bcf_read(test_bcf, test_header, test_record);

  while (chr_no < faidx_nseq(seq_ref)){
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    
    char *data = fai_fetch(seq_ref,name,&name_len);
    
    int vcf_pos = test_record->pos;
    std::cout << "VCF STARTT " << vcf_pos << std::endl;
    int start_pos = vcf_pos - 10;
    std::cout << "first START POS " << start_pos << std::endl;
    int end_pos = 19002055; //30001000
    char seqmod[1024];

    while(start_pos < end_pos){
      std::cout << "------ while ------" << std::endl;
      int readlength = drand48()*(80.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      std::cout << "start " << start_pos << " stop "<< stop << std::endl;
      
      if (start_pos <= vcf_pos)
      {
        std::cout << "VCF POS" << vcf_pos << std::endl;
        while (vcf_pos <= stop){
          bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL); 
          std::cout << "record"<< test_record->pos << std::endl;
          
          vcfstruct = bcf_read(test_bcf, test_header, test_record);
          vcf_pos = test_record->pos+1;
        }
      start_pos += readlength + 1;
      readlength = 0;
      }
    }
    /*while(bcf_read(test_bcf, test_header, test_record) == 0){
      std::cout << "Whilw "<< std::endl;
      int vcfpos = test_record->pos;
      bcf_read(test_bcf, test_header, test_record);
      vcfpos = test_record->pos;
      std::cout << "new pos "<< vcfpos << std::endl;
    }*/
  chr_no++;
  }
}

void VCF2(const char* fastafile,FILE *fp){  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;

  // initiate VCF
  htsFile *test_vcf = NULL;
  bcf_hdr_t *test_header = NULL;
  bcf1_t *test_record = bcf_init();
  //test_vcf = vcf_open("/home/wql443/WP1/data/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", "r");
  //test_vcf = vcf_open("/home/wql443/WP1/data/chr14_full.vcf", "r");
  test_vcf = bcf_open("/home/wql443/WP1/data/Chr14_subset.vcf", "r");
  test_header = bcf_hdr_read(test_vcf);

  while (chr_no < faidx_nseq(ref)){
    // The fasta file from the references
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);

    int start_pos =19000000; //20263521 ; //20250390;//20250390;// 19108144; 
    int end_pos = 19002055; //20266521 ;//;20255390;//20255390; //19111144;  

    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos_init = 0;
      int vcfpos = test_record->pos;

      while(start_pos < end_pos){
        std::srand(start_pos+std::time(nullptr));
        // Seed random number generator
        int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
        int frag_end = start_pos+rand_len;
        char* sequence = faidx_fetch_seq(ref,name,start_pos,frag_end,&name_len);
        char* pch;
        pch = strchr(sequence,'N');
        if (pch != NULL){
          //Disregards any read with 'N' in.. change this to just change the reading position
          start_pos += rand_len + 1;
        }
        else if (start_pos+1 <= vcfpos){
          std::cout << "------------ NEW REGION -------" << std::endl;          
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int length = strlen(sequence);
        
          //extend the sequence
          int alt_pos_ext=0;
          int ref_init_len=0;

          while (vcfpos <= frag_end){
            if (vcfpos == 19000655)
            {
              break;
            }
            else
            {
              if (test_record->n_allele > 1){
                std::cout << "---------- NEW ALTERNATIVE --------------" << std::endl;
                std::cout << "start " << start_pos << " vcf "<<  vcfpos << " end "<< frag_end << std::endl;
              }
              std::cout << Seq_str << std::endl;
              std::cout << " --- end of loop --- " << std::endl;
            }
            vcfpos_init = vcfpos;
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
            std::cout << vcfpos_init << " test2 " << vcfpos << std::endl;
            if(vcfpos_init == vcfpos){ start_pos += rand_len + 1; }
          }
        }
        else{
          std::cout << "else loop" << std::endl;
          std::cout << "init " << vcfpos_init << " pos " << vcfpos << std::endl;
          
          vcfpos_init = vcfpos;
            //iterating through the vcf file
          vcfstruct = bcf_read(test_vcf, test_header, test_record);
          vcfpos = test_record->pos+1;
          if(vcfpos_init == vcfpos){
            std::cout << "INIT LORT" << std::endl;
            break;}
        }
        start_pos += rand_len + 1;
      }
    }
    chr_no++;
  }
  bcf_hdr_destroy(test_header);
  bcf_destroy(test_record); 
  bcf_close(test_vcf);  
  return; 
}

static void delete_seq(char *str, size_t seq_len, size_t del_len, size_t pos)
{
    for (int i = 0; i < del_len; i++)
    {
        memmove(&str[pos], &str[pos+1], seq_len - pos);
        seq_len--;
    }
}

static void instert_seq(char *str, size_t len, char insert_seq[],size_t ins_len, size_t pos)
{
    for (int i = 0; i < ins_len; i++)
    {
        memmove(&str[pos+1], &str[pos], len - pos + 1);
        str[pos] = insert_seq[i];
        len++;
    }
    //the insertion doesn't overlap so i remove the nucleotide which the insertion replace
    delete_seq(str, strlen(str),1,pos+ins_len);
}

void variant_ref(const char* fastafile,FILE *fp){
  //creates fasta structure
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int chr_no = 0;
  
  // VCF STRUCTURE
  htsFile *test_bcf = NULL;
  bcf_hdr_t *test_header = NULL;
  //"/home/wql443/WP1/data/Chr14_subset.vcf"
  test_bcf = bcf_open("/home/wql443/WP1/data/test2.vcf", "r");
  test_header = bcf_hdr_read(test_bcf);
  bcf1_t *test_record = bcf_init();

  const char *name = faidx_iseq(seq_ref,chr_no);
  int name_len =  faidx_seq_len(seq_ref,name);
  //fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
  
  char *data = fai_fetch(seq_ref,name,&name_len);
  char seqmod[1024];

  kstring_t fa_kstr;
  fa_kstr.s = NULL; fa_kstr.l = fa_kstr.m = 0;

  /*std::cout << data[2] << std::endl;
  data[2] = 'X';
  std::cout << data[2] << std::endl;
  std::cout << strlen(data) << std::endl;*/

  /*14	19297364	.	TG	T
  14	19291757	.	G	GC
  14	19291758	.	AGTG	A
  14	19291763	.	T	C
  14	19291815	.	T	C
  14	19291845	.	T	C
  14	19291886	.	CTG	C
  14	19291899	.	C	CTGTT
  14	19325134	.	TA	TAA,T
  19020697	.	G	A,T
  14	19054940	.	A	AG,G
  */

  while(bcf_read(test_bcf, test_header, test_record) == 0){
    bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL);
    // bcf_hdr_id2name(test_header, test_record->rid) << "\t" <<
    int pos = (int) test_record->pos + 1;
    int n_allele = (int) test_record->n_allele;
    char* ref = test_record->d.allele[0];

    int readlength = drand48()*(80.0-30.0)+30.0;
    int idx = drand48()*(readlength); //where the variation is in the read
    
    int start = pos-idx-1;
    int stop = start + readlength;
    strncpy(seqmod,data+start,readlength);
    
    char seqbefore[1024];
    char seqafter[1024];
    char seq_indel[1024];
    //std::cout << n_allele << std::endl;

    if (n_allele < 3)
    {
      char* alt = test_record->d.allele[1];
      if (strlen(ref) == 1)
      {
        if (strlen(alt) == 1)
        {
          //std::cout << pos << "\t" << ref << "\t" << alt << std::endl;
          //std::cout << seqmod << std::endl;
          seqmod[idx] = *alt;    
          //std::cout << seqmod << std::endl;
          ksprintf(&fa_kstr,">%s:%d-%d_len:%d_%d_sub\n%s\n",name,start,stop,readlength,strlen(seqmod),seqmod);  
        }
        else //insertions
        {
          std::cout << pos << "\t" << ref << "\t" << alt << "\t index " << idx << std::endl;
          memcpy(seqbefore, &seqmod[0], idx); //from 1st position to position before the alternative start index
          memcpy(seqafter, &seqmod[idx+strlen(ref)], strlen(seqmod)); //from after the alternative to the end of the sequence
          std::cout << seqmod << std::endl;
          snprintf(seq_indel,1024,"%s%s%s\n",seqbefore,alt,seqafter);
          std::cout << seq_indel << std::endl;
          std::cout << "-------" << std::endl;
          
          //std::cout << seqmod << std::endl;
          //std::cout << seq_indel << std::endl;
          memset(seqbefore, 0, sizeof seqbefore);
          memset(seqafter, 0, sizeof seqafter);
          memset(seq_indel, 0, sizeof seq_indel);
        }
      }

      
    }
    
    /*
    if (pos == 21444573)
    {
      std::cout << pos << "\t" << ref << "\t" << alt << std::endl;
      std::cout << test_record->d.allele[1] << "\t" << test_record->d.allele[2] << std::endl;
      std::cout << n_allele << std::endl;
    }*/
    
    /*
    if (pos == 19291899)
    {
      int start = pos-idx;
      int stop = start + readlength;
      strncpy(seqmod,data+pos-idx,readlength);
      std::cout << seqmod << std::endl;
      std::cout << pos << "\t" << ref << "\t" << alt << std::endl;
      std::cout << data[pos-1] << std::endl;
      data[pos-1] = *alt;
      std::cout << data[pos-1] << std::endl;
      char seqbefore[1024];
      char seqafter[1024];
      memcpy(seqbefore, &seqmod[0], idx-1);
      memcpy(seqafter, &seqmod[idx+strlen(ref)-1], strlen(seqmod));
      std::cout << seqbefore << " " << seqafter << std::endl;
      char seq[1024];
      snprintf(seq,1024,"%s%s%s\n",seqbefore,alt,seqafter);
      std::cout << seq << std::endl;
      std::cout << " pos " << start << " end " << stop << " str len "<< strlen(seqmod) << strlen(seq) << std::endl;
      break;
    }*/
    memset(seqmod, 0, sizeof seqmod);
  }
  fwrite(fa_kstr.s,sizeof(char),fa_kstr.l,fp);fa_kstr.l =0;  
}

int main(int argc,char **argv){
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr14.fa";
  FILE *fp;
  fp = fopen("favcf.fa","wb");
  //VCF2(fastafile,fp); 
  srand48(time(NULL));
  variant_ref(fastafile,fp);
  fclose(fp);
  //std::cout << std::abs(10-5) << std::endl;

  return 0;
}
//g++ FaVCF_2.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl

/*while (chr_no < faidx_nseq(seq_ref)){
    //int whichref = lrand48() % faidx_nseq(seq_ref);
    //const char *name = faidx_iseq(seq_ref,whichref);
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);
    
    int start_pos = 20000000;
    int end_pos = 22000000;
    char seqmod[1024];

    while(start_pos <= end_pos){
      int readlength = drand48()*(80.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      
      //adress of entire array, first element always the same but adding start for query source start, and stop-start : number of elements
      strncpy(seqmod,data+start_pos,readlength);

      char * pch;
      pch = strchr(seqmod,'N');
      if (pch != NULL){start_pos += readlength + 1;}
      else {
        char nt[] = "tT";
        //Deamin_char(seqmod,nt,readlength);

        fprintf(fp,">%s:%d-%d_length:%d\n%s\n",name,start_pos,start_pos+readlength,readlength,seqmod);
        start_pos += readlength + 1;
        readlength = 0;
        memset(seqmod, 0, sizeof seqmod);
      }
    }
  chr_no++;
  }*/

/*
while(start_pos <= end_pos){
      int readlength = drand48()*(80.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      bcf_read(test_bcf, test_header, test_record);
      int vcf_cur =  test_record->pos;
      std::cout << "start " << start_pos << " end " << stop << std::endl;
      while(bcf_read(test_bcf, test_header, test_record) == 0) {
        std::cout << test_record->pos << std::endl;
      }
      /*while (start_pos < vcf_cur < stop)
      {
        std::cout << "vcfpos " << vcf_cur << std::endl;
        bcf_read(test_bcf, test_header, test_record);
        vcf_cur += test_record->pos;
        start_pos += readlength + 1;
      }
      
      / *bcf_read(test_bcf, test_header, test_record);
      std::cout << "vcfpos " << test_record->pos << std::endl;
      int vcfpos = test_record->pos;
      std::cout << "start " << start_pos << " end " << stop << std::endl;* /
      start_pos += readlength + 1;
      readlength = 0;* /
      start_pos += readlength + 1;
    }
*/