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

void VCF5(const char* fastafile,double cov=1.0){  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  std::ofstream outfa("hurra.fa");

  // initiate VCF
  htsFile *test_vcf = NULL;
  bcf_hdr_t *test_header = NULL;
  bcf1_t *test_record = bcf_init();
  //test_vcf = vcf_open("/home/wql443/WP1/data/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", "r");
  //test_vcf = vcf_open("/home/wql443/WP1/data/chr14_full.vcf", "r");
  test_vcf = vcf_open("/home/wql443/WP1/data/test2.vcf", "r");
  test_header = bcf_hdr_read(test_vcf);

  while (chr_no < faidx_nseq(ref)){
    // The fasta file from the references
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);

    int start_pos =1; //20263521 ; //20250390;//20250390;// 19108144; 
    int end_pos = name_len; //20266521 ;//;20255390;//20255390; //19111144;  

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
        int dist = init/cov * rand_len; 
        int frag_end = start_pos+rand_len;
        char* sequence = faidx_fetch_seq(ref,name,start_pos,frag_end,&name_len);
        char* pch;
        pch = strchr(sequence,'N');
        if (pch != NULL){
          //Disregards any read with 'N' in.. change this to just change the reading position
          start_pos += dist + 1;
        }
        else if (start_pos+1 <= vcfpos){
          std::cout << "------------ NEW REGION -------" << std::endl;
          char nt[] = "tT";
          
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int length = strlen(sequence);
        
          //extend the sequence
          int alt_pos_ext=0;
          int ref_init_len=0;

          while (vcfpos <= frag_end){
            bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL); // bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
            char* ref;
            char* alt;
            char* alt2;
            int alt_no;
            char CNV_char[] = "<";

            if (test_record->n_allele > 1){
              //NB! my alt_pos is an integer like 18, but it will be converted into an 0-index thus count 19 nt in the sequence
              int alt_pos;
              std::cout << "---------- NEW ALTERNATIVE --------------" << std::endl;
              std::cout << "start " << start_pos << " vcf "<<  vcfpos << " end "<< frag_end << std::endl;
              std::cout << Seq_str << std::endl;
              
              alt_pos = vcfpos-start_pos-1;
              if(vcfpos-start_pos == 0){std::cout << "LOORT" << std::endl;} //to ensure variation in the first basepair.. BUT THIS IS WRONG!!!

              ref = test_record->d.allele[0];
              alt = test_record->d.allele[1];

              alt_no = test_record->d.info->len;

              std::string refstr(ref);
              std::string altstr(alt);

              int refsize = (int) refstr.size();
              int altsize = (int) altstr.size();
              int new_alt_pos = alt_pos+alt_pos_ext;
              if (alt_no < 2){
                std::cout << "single alternative ref " << ref << " alt " << altstr << std::endl;
                //single alternative
                ref_init_len = refsize;
                
                if (new_alt_pos<0){
                  //deletions from previous alternative have removed this current alternative
                  Seq_str = Seq_str;
                  std::cout << "only if alt pos" << new_alt_pos << std::endl;
                  //Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                }
                //variation within the "middle"
                else if (0<=new_alt_pos && new_alt_pos<Seq_str_alt.length()){
                  std::cout << "else if alt pos" << new_alt_pos << std::endl;
                  Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                }
                //No variation
                else
                {
                  std::cout << "only else  alt pos" << new_alt_pos << std::endl;
                  Seq_str = Seq_str;
                }
                //update the index to consider the current variation changes
                if(refsize == altsize){
                  //substituion - 0 change
                  alt_pos_ext += (refsize-altsize);
                }
                else if (refsize > altsize){
                  //deletion - index are decreased
                  alt_pos_ext -= (refsize-altsize);
                }
                else if (altsize > refsize){
                  //insertions - index are increased
                  alt_pos_ext += std::abs(refsize-altsize);
                }
              }
              else{
                // insert the removed parts here for multiple alleles
                vcfstruct = bcf_read(test_vcf, test_header, test_record);
                vcfpos = test_record->pos+1;
              }
            std::cout << Seq_str << std::endl;
            std::cout << " --- end of loop --- " << std::endl;
            }
            vcfpos_init = vcfpos;
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
            std::cout << vcfpos_init << " test2 " << vcfpos << std::endl;
            if(vcfpos_init == vcfpos){ start_pos += dist + 1; }
          }

          outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
          outfa << Seq_str << std::endl;  
        }
        else{
          std::cout << "else loop" << std::endl;
          while (vcfpos < start_pos){
            vcfpos_init = vcfpos;
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
            if(vcfpos_init == vcfpos){ start_pos += dist + 1;; }
          }
        }
        start_pos += dist + 1;
      }
    }
    chr_no++;
  }
  bcf_hdr_destroy(test_header);
  bcf_destroy(test_record); 
  bcf_close(test_vcf);  
  return; 
}

int main(int argc,char **argv){
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr14.fa";
  VCF5(fastafile);
  //std::cout << std::abs(10-5) << std::endl;
  return 0;
}

//g++ FaVCF_thread.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl


/*
//selecting one random of the multiple random alleles
                int rand_alt = rand() % alt_no + 1;
                std::cout << "RANDOM" << rand_alt << " alt " << alt_no << std::endl;

                alt2 = test_record->d.allele[rand_alt];
                std::string altstr2(alt2);
                int alt2size = (int) altstr2.size();  
                std::cout << "mulitple alternatives "<< std::endl;
                std::cout << "seq length" << Seq_str.length() << std::endl;
                std::cout << "ref " << ref << " size " << refsize << std::endl;
                std::cout << "alt1 " << alt << " alt2 "<< alt2 << std::endl;
                
                if (alt[0] == '<'){
                  //CNV variation starts with '<'
                  altstr2.erase(0,3).pop_back();
                  int CNV_2 = stoi(altstr2);
                  
                  //creates the copy number string
                  std::string CNV_2_str(CNV_2, ref[0]); 

                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),CNV_2_str);
                }
                else
                {  
                  // if the alternative alleles are on the first basepair the index are affected as its normally both insertion and deltion of the ref
                  if (vcfpos-start_pos == 0) //keep index in seq for variation if first basepar
                  { 
                    std::cout <<  "multiple if "<< std::endl;    
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else if ((vcfpos-refsize)-start_pos < 0) //variation in the first basepairs with one alternative -> deletions removing index 0
                  {
                    std::cout <<  "multiple else if "<< std::endl;
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else
                  {
                    std::cout <<  "multiple else "<< std::endl;
                    //multiple different unique alternatives
                    Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr2);
                  }
                  if(refsize == alt2size){alt_pos_ext += 0;}
                  else if (refsize > alt2size){alt_pos_ext -= (refsize-alt2size);}
                  else if (alt2size > refsize){alt_pos_ext += (alt2size-refsize);}
                }*/