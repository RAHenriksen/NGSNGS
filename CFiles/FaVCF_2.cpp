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
  // "/home/wql443/WP1/data/Chr14_subset.vcf  "/home/wql443/WP1/data/test2.vcf"   "/home/wql443/WP1/data/chr14_full.vcf"
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

  while(bcf_read(test_bcf, test_header, test_record) == 0){
    bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL);
    // bcf_hdr_id2name(test_header, test_record->rid) << "\t" <<
    int pos = (int) test_record->pos + 1;
    int n_allele = (int) test_record->n_allele;
    char* ref = test_record->d.allele[0];
    char* alt = test_record->d.allele[1];

    int readlength = drand48()*(80.0-30.0)+30.0;
    int idx = drand48()*(readlength); //where the variation is in the read
    
    int start = pos-idx-1;
    int stop = start + readlength;
    strncpy(seqmod,data+start,readlength);
    
    char seqbefore[1024];
    char seqafter[1024];
    char seq_indel[1024];
    //std::cout << n_allele << std::endl;

    int rand_alt = rand() % (n_allele-1) +1;
    char* alt_select = test_record->d.allele[rand_alt];

    if (alt[0] != '<'){ // canonical variations
      if (strlen(ref) == strlen(alt_select)){
        memcpy(seqbefore, &seqmod[0], idx); //from 1st position to position before the alternative start index
        memcpy(seqafter, &seqmod[idx+strlen(ref)], strlen(seqmod)); //from after the alternative to the end of the sequence
        snprintf(seq_indel,1024,"%s%s%s\n",seqbefore,alt_select,seqafter);
        ksprintf(&fa_kstr,">%s:%d-%d_%d_len:%d_%d_sub\n%s\n",name,start+1,stop,start+idx+1,strlen(seqmod),strlen(seq_indel)-1,seq_indel);
      }
      else{ //indels
        memcpy(seqbefore, &seqmod[0], idx); //from 1st position to position before the alternative start index
        memcpy(seqafter, &seqmod[idx+strlen(ref)], strlen(seqmod)); //from after the alternative to the end of the sequence
        snprintf(seq_indel,1024,"%s%s%s\n",seqbefore,alt_select,seqafter);
        if ((strlen(seq_indel)-1) > strlen(seqmod)){
          //insertions
          ksprintf(&fa_kstr,">%s:%d-%d_%d_len:%d_%d_ins\n%s\n",name,start+1,stop,start+idx+1,strlen(seqmod),strlen(seq_indel)-1,seq_indel);
        }
        else{
          //deletions
          ksprintf(&fa_kstr,">%s:%d-%d_%d_len:%d_%d_del\n%s\n",name,start+1,stop,start+idx+1,strlen(seqmod),strlen(seq_indel)-1,seq_indel);
        } 
      }
    }
    else{ //Copy Number Variations
      std::cout << "else loop" << std::endl;
      int cnv_val = atoi(&test_record->d.allele[rand_alt][3]);
      std::cout << "random val" << rand_alt << std::endl;
      char* alt_select = test_record->d.allele[rand_alt];
      std::cout << "random char " << alt_select << std::endl;
      std::cout << cnv_val << std::endl;
      char cnv_char[1024];
      std::cout << "ref 1 " << ref << std::endl;
      memset(cnv_char, *ref, cnv_val);
      std::cout << "ref 2 " << cnv_char << std::endl;
      cnv_char[cnv_val] = '\0';

      memcpy(seqbefore, &seqmod[0], idx); //from 1st position to position before the alternative start index
      memcpy(seqafter, &seqmod[idx+strlen(ref)], strlen(seqmod)); //from after the alternative to the end of the sequence
      snprintf(seq_indel,1024,"%s%s%s\n",seqbefore,cnv_char,seqafter);

      if ((strlen(seq_indel)-1) > strlen(seqmod)){
        //insertions
        ksprintf(&fa_kstr,">%s:%d-%d_%d_len:%d_%d_ins\n%s\n",name,start+1,stop,start+idx+1,strlen(seqmod),strlen(seq_indel)-1,seq_indel);
      }
      else{
        //deletions
        ksprintf(&fa_kstr,">%s:%d-%d_%d_len:%d_%d_del\n%s\n",name,start+1,stop,start+idx+1,strlen(seqmod),strlen(seq_indel)-1,seq_indel);
      } 
    }

    memset(seqmod, 0, sizeof seqmod);
    memset(seqbefore, 0, sizeof seqbefore);
    memset(seqafter, 0, sizeof seqafter);
    memset(seq_indel, 0, sizeof seq_indel);
  }
  fwrite(fa_kstr.s,sizeof(char),fa_kstr.l,fp);fa_kstr.l =0;  
}

int main(int argc,char **argv){
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr14.fa";
  FILE *fp;
  fp = fopen("favcf.fa","wb");
  //VCF2(fastafile,fp); 
  srand48(time(NULL));
  srand(time(NULL));
  variant_ref(fastafile,fp);
  fclose(fp);
  //std::cout << std::abs(10-5) << std::endl;

  return 0;
}
//g++ FaVCF_2.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
