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

char* HaploGenome(char* genome,char genome_data1[],char genome_data2[],int chr_sizes,const char* bcf_file,const char *chr_names[],const char* VarType,FILE *VarInfoFile,const char* HeaderIndiv){

  DNA_CAPITAL(genome_data1);DNA_CAPITAL(genome_data2);

  htsFile *bcf_obj = NULL;
  bcf_hdr_t *bcf_head = NULL;
  // "/home/wql443/WP1/data/Chr14_subset.vcf  "/home/wql443/WP1/data/test2.vcf"   "/home/wql443/WP1/data/chr14_full.vcf"
  bcf_obj = bcf_open(bcf_file, "r");
  hts_idx_t *idx = bcf_index_load(bcf_file);
  bcf_head = bcf_hdr_read(bcf_obj);
  fprintf(stderr, "File contains before %i samples\n",bcf_hdr_nsamples(bcf_head)); //bcf.nsamples(fp)
  int nseq = 0;

  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(bcf_head, &nseq);

  // check if the headers match
  int chr_exit = 0;
  for (int vcf_i = 0; vcf_i < nseq; vcf_i++){
    if(strcasecmp(seqnames[vcf_i],chr_names[0])==0){continue;}//fprintf(stderr,"number of chromosomes %d and names %s\n",chr_i,chr_names[chr_i]);
    else{chr_exit++;}
  }
  if(chr_exit == nseq){
    fprintf(stderr,"Reference file chromosome %s is not identified in vcf/bcf header\n",chr_names[0]);
    exit(0);
  }
  
  // Set sample number 
  //const char* HeaderIndiv = "HG00096";
  
  //bcf_hdr_t *bcf_head2 =bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);
  
  //int indiv = bcf_hdr_set_samples(bcf_head, "HG00099", 0);
  //fprintf(stderr, "File contains after %i samples %d\n",bcf_hdr_nsamples(bcf_head),indiv); //bcf.nsamples(fp)
  //bcf_hdr_parse_sample_line(bcf_head,htxt)
  //initializing variation sampling 
  bcf1_t *bcf_records = bcf_init();
  if (bcf_records == NULL){fprintf(stderr,"WARNING NO VARIANTS IN VCF FILE"); exit(0);}
  else if(bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0){
    // else if{bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0; 
    //bcf_read(bcf_obj, bcf_head, bcf_records)==0; //ignoring return value (int)
    
    /*kstring_t shared,indiv = bcf_records->indiv;
    std::cout << shared << std::endl;
    fprintf(stderr,"test %s \n",indiv.s);
    exit(0);*/

    /*
    For indels, which has previously been removed and not yet incorporated and improved
    size_t chr_pos_before;
    size_t chr_pos_after;
    size_t insert_total = 0;
    size_t del_total = 0;
    size_t old_pos = 0;
    */

    hts_itr_t *itr = bcf_itr_querys(idx, bcf_head,chr_names[0]);
    int record_indiv;
    while ((record_indiv = bcf_itr_next(bcf_obj, itr, bcf_records)) == 0){
      //fprintf(stderr,"WHILE LOOP\n");
      if(strcasecmp(bcf_hdr_id2name(bcf_head, bcf_records->rid),chr_names[0])==0){
        bcf_unpack((bcf1_t*)bcf_records, BCF_UN_ALL);

        //https://github.com/samtools/htslib/blob/226c1a813bc5d0582f7e0b0bdb4b3ea9e3ee4ce4/htslib/vcf.h#L1019
        int nsamples = bcf_hdr_nsamples(bcf_head);

        //fprintf(stderr,"number of samples %d\n",nsamples);

        // gt data for each call
        int32_t ngt_arr = 0;     
        int32_t *gt_arr = NULL;
        int ngt = bcf_get_genotypes(bcf_head, bcf_records, &gt_arr, &ngt_arr);
        int max_ploidy = ngt/nsamples;

        
        for (int i =0; i<nsamples; i++){
          char* haplotype1 = NULL;char* haplotype2 = NULL;//char haplotype1[1024] = {0};char haplotype2[1024] = {0}; //
          //iterates through all samples
          // fprintf(stderr,"--------\nTHE SAMPLE INDEX IS %d and sample name is %s\n",i,bcf_head->samples[i]);
          // match the names in the header with the input indivduals
          if(strcasecmp(bcf_head->samples[i],HeaderIndiv)==0){
            // fprintf(stderr,"--------\nTHE SAMPLE INDEX IS %d and sample name is %s and input %s\n",i,bcf_head->samples[i],HeaderIndiv);
            int32_t *ptr = gt_arr + i*max_ploidy;
            int haplotype1_int; int haplotype2_int; 
            //for (int j=0; j<max_ploidy; j++){fprintf(stderr,"the ploidy index is %d \t alelle values %d\n",j,bcf_gt_allele(ptr[j]));}
            haplotype1_int = bcf_gt_allele(ptr[0]); haplotype2_int = bcf_gt_allele(ptr[1]);

            //Extract actual allele chars;
            if(haplotype1_int == 0 && haplotype2_int == 0){
              //homozygous for reference
              //fprintf(stderr,"homozygous for reference\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
            }
            else if(haplotype1_int == 1 && haplotype2_int == 1){
              //homozygous for alternative
              //fprintf(stderr,"homozygous for alternative\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
            }
            else if(haplotype1_int != haplotype2_int){
              //heterozygous
              //fprintf(stderr,"heterozygous\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
            }

            size_t pos = (int) bcf_records->pos;
            //fprintf(stderr,"Position %zu \t The haplotypes values are %d \t %d and characters are %c \t %c \n",pos,haplotype1_int,haplotype2_int,*haplotype1,*haplotype2);
            //Current SNP location
            //fprintf(stderr,"before alterations %c%c%c\t%c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1],genome_data2[pos-1],genome_data2[pos],genome_data2[pos+1]);
            genome_data1[pos] = *haplotype1;genome_data2[pos] = *haplotype2;
            //fprintf(stderr,"after alterations %c%c%c\t%c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1],genome_data2[pos-1],genome_data2[pos],genome_data2[pos+1]);
            
          }
        }
      }
      else{
        fprintf(stderr,"Reference chromosome %s has no variations in vcf file at chromosomes %s\n",chr_names[0],bcf_hdr_id2name(bcf_head, bcf_records->rid));
      }
    }
  }
  //valgrind ./ngsngs -i ../vcfdata/chr14.fa -r 10000000 -t1 4 -l 100 -seq SE -f fq -ne -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -bcf ../vcfdata/chr14_indiv_3.bcf -v snp -indiv HG00099 -o Chr14_HG00099
  bcf_hdr_destroy(bcf_head);
  bcf_destroy(bcf_records); 
  bcf_close(bcf_obj);
  //fprintf(stderr,"Genome 1    example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  sprintf(genome+strlen(genome_data1),"%s",genome_data1);
  strcat(genome,genome_data2);
  /*fprintf(stderr,"Haplotype 1 example %c%c%c%c%c\n",genome_data1[19000014],genome_data1[19000015],genome_data1[19000016], genome_data1[19000017], genome_data1[19000018]);
  fprintf(stderr,"Haplotype 2 example %c%c%c%c%c\n",genome_data2[19000014],genome_data2[19000015],genome_data2[19000016], genome_data2[19000017], genome_data2[19000018]);
  fprintf(stderr,"Genome II   example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  fprintf(stderr,"Genome III  example %c%c%c%c%c\n",genome[19000014+107349540],genome[19000015+107349540],genome[19000016+107349540], genome[19000017+107349540], genome[19000018+107349540]);*/

  return genome;
}



char* HaploType(char genome_data1[],htsFile *bcf_obj,bcf_hdr_t *bcf_head,hts_idx_t *bcf_idx,bcf1_t *bcf_records,const char *chr_name,int HaploGroup,const char* HeaderIndiv){

  fprintf(stderr, "INSIDE HAPLOTYPE FUNC \n"); //bcf.nsamples(fp)
  DNA_CAPITAL(genome_data1);
  std::cout << genome_data1[10000000] << std::endl;
  
  hts_itr_t *itr = bcf_itr_querys(bcf_idx, bcf_head,chr_name);
  int record_indiv;
  int nsamples = bcf_hdr_nsamples(bcf_head);
  fprintf(stderr,"number of samples %d\n",nsamples);
  
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
          //Current SNP location
          //fprintf(stderr,"%zu POSITION \t before alterations %c%c%c\n",pos,genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1]);
          genome_data1[pos] = *haplotype1;
          //fprintf(stderr,"after alterations %c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1]);
        }
      }
    }
    else{
      fprintf(stderr,"LOL\n");
      //fprintf(stderr,"Reference chromosome %s has no variations in vcf file at chromosomes %s\n",chr_names[0],bcf_hdr_id2name(bcf_head, bcf_records->rid));
    }
  }
  return genome_data1;
}



char* full_vcf_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[],const char* bcf_file,const char* VarType,const char* HeaderIndiv){
  fprintf(stderr,"BCF FILE %s\n",bcf_file);
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  
  const char *chr_name = faidx_iseq(seq_ref,0);
  int chr_len = faidx_seq_len(seq_ref,chr_name);
  chr_sizes[0] = chr_len;chr_sizes[1] = chr_len;
  chr_names[0] = chr_name;chr_names[1] = chr_name;
  genome_size += chr_len;
  chr_size_cumm[1] = genome_size;
  genome_size += chr_len;
  chr_size_cumm[2] = genome_size;

  char* genome = (char*) malloc(sizeof(char) * (genome_size+(chr_total*2))*2);
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  FILE *VarInfoFile = fopen("Variant_pos_log.txt", "w");
  fprintf(VarInfoFile,"Type\tChr\tRef_Pos\tRef\tRead_Pos\tAlt\n");
  for (int i = 0; i < chr_total; i++){

    char *Chr_hapl1 = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    char *Chr_hapl2 = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);

    if (Chr_hapl1 != NULL){
      fprintf(stderr,"INSIDE NULL DATA \t WITH VARTYPE %s and chr %s\n",VarType,chr_names[i]);
      //char* HaploGenome(char genome_data[],int chr_sizes,const char* bcf_file,const char *chr_names[]){
      
      HaploGenome(genome,Chr_hapl1, Chr_hapl2,chr_sizes[i],bcf_file,&chr_names[i],VarType,VarInfoFile,HeaderIndiv); //Chr14_subset4_bcf
    }
    //std::cout<< genome[19000016] << genome[19000017] << genome[19000018] << genome[19000016+107349540] << genome[19000017+107349540] << genome[19000018+107349540] << std::endl;
    //std::cout << strlen(genome)<<std::endl;
    // several of the build in functions allocates memory without freeing it again.
    free((char*)Chr_hapl1); //Free works on const pointers, so we have to cast into a const char pointer
    free((char*)Chr_hapl2); //Free works on const pointers, so we have to cast into a const char pointer
  }
  /*fprintf(stderr,"Genome IIII example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  fprintf(stderr,"Genome IV   example %c%c%c%c%c\n",genome[19000014+107349540],genome[19000015+107349540],genome[19000016+107349540], genome[19000017+107349540], genome[19000018+107349540]);
  fprintf(stderr,"DONE WITH DATA FOR LOOP \n");*/
  fclose(VarInfoFile);
  return genome;
}

