#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::reverse
#include <iostream>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/cram.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>
#include <random>
#include <map>
#include <math.h>
#include <pthread.h>

#include "NGSNGS_func.h"
#include "NtSubModels.h"
#include "mrand.h"
#include "RandSampling.h"

//#if defined(__APPLE__) && defined(__MACH__) 
//#include "NGSNGS_Random.h"
//#endif /* __APPLE__ */

#define LENS 4096
#define MAXBINS 100

void DNA_CAPITAL(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'A';
        break;
      case 'G':
      case 'g':
        *seq = 'G';
        break;
      case 'C':
      case 'c':
        *seq = 'C';
        break;
      case 'T':
      case 't':
        *seq = 'T';
        break;
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
  }
}

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void ReversComplement(char seq[]){
  // generates the reverse complementary sequence from an input sequence
  unsigned char nuc2int[255];
  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4;

  char ntdeam[4] = {'T', 'A', 'C', 'G'};
  char seq_intermediate[1024] = {0};
  strcpy(seq_intermediate,seq);
  //fprintf(stderr,"SEQUENCE \t\t%s\n",seq_intermediate);
  int seqlen = strlen(seq);
  //Complementing sequence
  for(int i=0;i<seqlen;i++){
    seq_intermediate[i] = ntdeam[nuc2int[(unsigned char) seq_intermediate[i]]]; //warning: array subscript has type 'char' [-Wchar-subscripts]
  }
  //fprintf(stderr,"COMP SEQUENCE \t\t%s\n",seq_intermediate);

  //reverse complement
  for(int i=seqlen-1;i>-1;i--){
    seq[seqlen-i-1] = seq_intermediate[i];
  }
  //just to ensure no issues arise in case of not clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,fasta_sampler *fs,char CommandArray[1024],const char* version){
  // Creates a header for the bamfile. The header is initialized before the function is called //
  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}

  // Creating header information
  
  char genome_len_buf[1024];
  //sam_hdr_add_line(header, "HD", "VN",version, "SO", "unsorted", NULL);
  for(int i=0;i<fs->nref;i++){
    snprintf(genome_len_buf,1024,"%d", fs->seqs_l[i]);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", fs->seqs_names[i], "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
   
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  
  }
  // Adding PG tag specifying the command used for simulations
  sam_hdr_add_pg(header,"NGSNGS","VN",version,"CL",CommandArray,NULL);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

char* HaploType(char genome_data1[],htsFile *bcf_obj,bcf_hdr_t *bcf_head,hts_idx_t *bcf_idx,bcf1_t *bcf_records,const char *chr_name,int HaploGroup,const char* HeaderIndiv){
  fprintf(stderr, "INSIDE HAPLOTYPE FUNC \n"); //bcf.nsamples(fp)
  DNA_CAPITAL(genome_data1);
  
  hts_itr_t *itr = bcf_itr_querys(bcf_idx, bcf_head,chr_name);
  int record_indiv;
  int nsamples = bcf_hdr_nsamples(bcf_head);
  
  while ((record_indiv = bcf_itr_next(bcf_obj, itr, bcf_records)) == 0){
    if(strcasecmp(bcf_hdr_id2name(bcf_head, bcf_records->rid),chr_name)==0){
      bcf_unpack((bcf1_t*)bcf_records, BCF_UN_ALL);

      // gt data for each call
      int32_t ngt_arr = 0;     
      int32_t *gt_arr = NULL;
      int ngt = bcf_get_genotypes(bcf_head, bcf_records, &gt_arr, &ngt_arr);
      int max_ploidy = ngt/nsamples;
      int HaploGroupAlter = max_ploidy-HaploGroup-1;

      //fprintf(stderr,"HaploGroup 1 %d and 2 %d\n",HaploGroup,HaploGroupAlter);

      for (int i =0; i<nsamples; i++){
        char* haplotype1 = NULL;
        //iterates through all samples
        // match the names in the header with the input indivduals
        if(strcasecmp(bcf_head->samples[i],HeaderIndiv)==0){
          //fprintf(stderr,"--------\nTHE SAMPLE INDEX IS %d and sample name is %s\n",i,bcf_head->samples[i]);
          int32_t *ptr = gt_arr + i*max_ploidy;

          //for (int j=0; j<max_ploidy; j++){fprintf(stderr,"the ploidy index is %d \t alelle values %d\n",j,bcf_gt_allele(ptr[j]));}
          int haplotype1_int = bcf_gt_allele(ptr[HaploGroup]); 
          int haplotype2_int = bcf_gt_allele(ptr[HaploGroupAlter]);
          
          //Extract actual allele chars;
          if(haplotype1_int == 0 && haplotype2_int == 0){
            //homozygous for reference
            haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
          }
          else if(haplotype1_int == 1 && haplotype2_int == 1){
            //homozygous for alternative
            haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
          }
          else if(haplotype1_int != haplotype2_int){
            //heterozygous
            haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[HaploGroup])];
          }

          size_t pos = (int) bcf_records->pos;
          //fprintf(stderr,"Position %zu \t The haplotypes values are %d \t %d and characters are %c \t %c \n",pos,haplotype1_int,haplotype2_int,*haplotype1,*haplotype2);
          //fprintf(stderr,"%zu POSITION \t before alterations %c%c%c\n",pos,genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1]);
          genome_data1[pos] = *haplotype1;
          //fprintf(stderr,"after alterations %c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1]);
        }
      }
    }
    else{
      fprintf(stderr,"Reference chromosome %s has no variations in vcf file at chromosomes %s\n",chr_name,bcf_hdr_id2name(bcf_head, bcf_records->rid));
    }
  }
  return genome_data1;
}