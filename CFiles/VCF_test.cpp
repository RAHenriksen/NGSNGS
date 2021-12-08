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

void VCF(const char* fastafile,double cov=1.0){  

  // Tried to use vector to store the different version of the alternatives.

  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  std::ofstream outfa("test.fa");

  // initiate VCF
  htsFile *test_vcf = NULL;
  bcf_hdr_t *test_header = NULL;
  bcf1_t *test_record = bcf_init();
  //test_vcf = vcf_open("/home/wql443/WP1/data/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", "r");
  test_vcf = vcf_open("/home/wql443/WP1/data/chr14_full.vcf", "r");
  test_header = bcf_hdr_read(test_vcf);

  while (chr_no < faidx_nseq(ref)){
    // The fasta file from the references
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);

    int start_pos = 19460500;  // 19060000;// 19135000; //19000000; //19066080;
    int end_pos = 19461000; // 19067000;//19136000; //19001000; //19088300;

    //14	19135717	.	GA	AA,G	100 //19135017; //19136017
    //14	19237205	.	T	TATGTTATG	100  //19236505;  //19238005;
    //14	19291899	.	C	CTGTT	100 //19291650; //19291900; -> CTGTTG


    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    //std::cout << "ts" << chrno << std::endl;
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos = test_record->pos+1;
      //std::cout << "vcf pos " << vcfpos << std::endl; 

      while(start_pos < end_pos){
        std::srand(start_pos+std::time(nullptr));
        // Seed random number generator
        int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
        //std::cout << " random " << rand_len << std::endl;
        int dist = init/cov * rand_len; 
        
        char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
        char* pch;
        pch = strchr(sequence,'N');
        if (pch != NULL){
          //Disregards any read with 'N' in.. change this to just change the reading position
          start_pos += dist + 1;
        }
        else if (start_pos <= vcfpos){
          // i'm not putting the deamination inside the while loop since that would for each ALT perform deamination
          //char nt[] = "tT";
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int alt_pos_etd=0;
          int alt2_pos_etd=0;
          
          //std::cout << "chr pos " << start_pos << " vcf pos " << vcfpos << std::endl;
          int length = strlen(sequence);
          
          std::vector<std::string> Seq_vec;
          Seq_vec.resize(30);
          Seq_vec.push_back(Seq_str);
          
          while (vcfpos <= start_pos+rand_len){
            bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL); // bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
            char* ref;
            char* alt;
            char* alt2;

            if (test_record->n_allele > 1){
              //NB! my alt_pos is an integer like 18, but it will be converted into an 0-index thus count 19 nt in the sequence
              int alt_pos;
              
              if (vcfpos-start_pos > 0){
                alt_pos = vcfpos-start_pos-1;
              }
              else
              {
                alt_pos = vcfpos-start_pos;
              }

              ref = test_record->d.allele[0];
              alt = test_record->d.allele[1];
              alt2 = test_record->d.allele[2];

              std::string refstr(ref);
              std::string altstr(alt);
              int refsize = (int) refstr.size();
              int altsize = (int) altstr.size();
              
              if (alt2==NULL){
                //single alternative allele
                for (int i = 0; i < Seq_vec.size(); i++){Seq_vec[i] = Seq_vec[i].replace(alt_pos+alt_pos_etd,strlen(ref),altstr);}
              }
              else if(alt2!=NULL&&(strlen(alt2)==0)){
                for (int i = 0; i < Seq_vec.size(); i++){Seq_vec[i] = Seq_vec[i].replace(alt_pos+alt_pos_etd,strlen(ref),altstr);}
              }
              else if(alt2!=NULL&&(strlen(alt2)>0)){
                std::string altstr2(alt2);
                int i=2;
                std::vector<std::string> Alt_tmp;
                //Seq_vec.push_back(Seq_str);
                //Seq_vec[0].replace(alt_pos+alt2_pos_etd,strlen(ref),altstr2);
                while (test_record->d.allele[i] != NULL){
                  //std::string test(test_record->d.allele[2]);
                  //std::cout << "test " << test << std::endl;
                  for (int j = 0; j < Seq_vec.size(); j++){
                    std::string Seq_tmp = Seq_vec[j];
                    //std::cout << "orig "<< Seq_str << std::endl;
                    //std::cout << "temp " << Seq_tmp << std::endl;
                    Alt_tmp.push_back(Seq_tmp.replace(alt_pos+alt2_pos_etd,strlen(ref),(std::string) test_record->d.allele[i]));
                  }
                  i++;
                  //std::cout<< (std::string) test_record->d.allele[i] << std::endl;
                  //std::cout << "size " <<Alt_tmp.size() << std::endl;
                  //for (int i = 0; i < Alt_tmp.size(); i++){std::cout << "tmp "<<Alt_tmp[i] << "\n";}
                }
                Seq_vec = Alt_tmp;
                
                //multiple alternative
                
                int alt2size = (int) altstr2.size();
                if(refsize > alt2size){alt2_pos_etd = alt2_pos_etd - std::abs(refsize-alt2size);}
                else if(refsize <= alt2size){alt2_pos_etd = alt2_pos_etd + std::abs(refsize-alt2size);}
                //alt2_pos_etd
                
                //std::cout << "multiple alternative" << std::endl;
                //std::cout << "pos " << alt_pos << " ref " << refstr << " alt " << altstr << " alt2 " << altstr2 << std::endl;

                //Seq_str = Seq_str.replace(alt_pos+alt_pos_etd,strlen(ref),altstr);
                //Seq_str_alt = Seq_str_alt.replace(alt_pos+alt2_pos_etd,strlen(ref),altstr2);
                //std::cout << "---" << std::endl;
              }
            
              if(refsize > altsize){alt_pos_etd = alt_pos_etd - std::abs(refsize-altsize);}
              else if(refsize <= altsize){alt_pos_etd += std::abs(refsize-altsize);}

              //Conver char* to strings in order to perform SNV or indels based on the reference and alternative.

            }
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }

          for (int i = 0; i < Seq_vec.size(); i++){
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_vec[i].length() << std::endl;
            outfa << Seq_vec[i] << std::endl;
          }
        }
        else{
          while (vcfpos < start_pos){
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }
        }
        start_pos += dist + 1;
      }
    }
    chr_no++;
  }
  //bcf_hdr_destroy(test_header);
  //bcf_destroy(test_record); 
  //bcf_close(test_vcf);  
  return; 
}


void VCF2(const char* fastafile,double cov=1.0){  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  std::ofstream outfa("test.fa");

  // initiate VCF
  htsFile *test_vcf = NULL;
  bcf_hdr_t *test_header = NULL;
  bcf1_t *test_record = bcf_init();
  //test_vcf = vcf_open("/home/wql443/WP1/data/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz", "r");
  test_vcf = vcf_open("/home/wql443/WP1/data/chr14_full.vcf", "r");
  test_header = bcf_hdr_read(test_vcf);

  while (chr_no < faidx_nseq(ref)){
    // The fasta file from the references
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);

    int start_pos = 1;  // 19060000;// 19135000; //19000000; //19066080;
    int end_pos = name_len; // 19067000;//19136000; //19001000; //19088300;

    //14	19135717	.	GA	AA,G	100 //19135017; //19136017
    //14	19237205	.	T	TATGTTATG	100  //19236505;  //19238005;
    //14	19291899	.	C	CTGTT	100 //19291650; //19291900; -> CTGTTG

    //alternative number 5 vcf pos 34215377     TTTATTTA TTTTATTTATTTATTTATTTA T

    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    //std::cout << "ts" << chrno << std::endl;
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos = test_record->pos+1;
      //std::cout << "vcf pos " << vcfpos << std::endl; 

      while(start_pos < end_pos){
        std::srand(start_pos+std::time(nullptr));
        // Seed random number generator
        int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
        //std::cout << " random " << rand_len << std::endl;
        int dist = init/cov * rand_len; 
        
        char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
        char* pch;
        pch = strchr(sequence,'N');
        if (pch != NULL){
          //Disregards any read with 'N' in.. change this to just change the reading position
          start_pos += dist + 1;
        }
        else if (start_pos <= vcfpos){
          // i'm not putting the deamination inside the while loop since that would for each ALT perform deamination
          //char nt[] = "tT";
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int alt_pos_etd=0;
          int alt2_pos_etd=0;
          
          //std::cout << "chr pos " << start_pos << " vcf pos " << vcfpos << std::endl;
          int length = strlen(sequence);

          while (vcfpos <= start_pos+rand_len){
            bcf_unpack((bcf1_t*)test_record, BCF_UN_ALL); // bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
            char* ref;
            char* alt_1;
            int alt_no;

            if (test_record->n_allele > 1){
              //NB! my alt_pos is an integer like 18, but it will be converted into an 0-index thus count 19 nt in the sequence
              int alt_pos;
              
              if (vcfpos-start_pos > 0){
                alt_pos = vcfpos-start_pos-1;
              }
              else
              {
                alt_pos = vcfpos-start_pos;
              }

              ref = test_record->d.allele[0];
              alt_1 = test_record->d.allele[1];
              alt_no = test_record->d.info->len;

              char test[] = "<";

              char* alt_2;
              alt_2 = test_record->d.allele[alt_no];
              std::string alt_no1(alt_1);
              std::string alt_no2(alt_2);

              if (alt_no > 1){
                std::cout << "alternative number " << alt_no << " vcf pos "<< vcfpos << std::endl;
                if (alt_1[0] == '<'){
                  //removes <CN from all the copy number alternatives
                  alt_no1.erase(0,3).pop_back();  
                  alt_no2.erase(0,3).pop_back();
                  int CNV_1 = stoi(alt_no1);
                  int CNV_2 = stoi(alt_no2);
                  
                  std::string CNV_1_str(CNV_1, ref[0]); 
                  std::string CNV_2_str(CNV_2, ref[0]); 

                  //std::cout << Seq_str << std::endl;
                  Seq_str = Seq_str.replace(alt_pos+alt_pos_etd,strlen(ref),CNV_1_str);
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_etd,strlen(ref),CNV_2_str);
                  //std::cout << Seq_str << std::endl;
                  //std::cout << Seq_str_alt << std::endl;
                  //std::cout << test_record->d.allele[1] << " " << test_record->d.allele[4][1] << " " << test_record->d.allele[5] << std::endl;
                }
                else
                {
                  Seq_str = Seq_str.replace(alt_pos+alt_pos_etd,strlen(ref),alt_no1);
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_etd,strlen(ref),alt_no2);
                  //Seq_str = Seq_str.replace(alt_pos+alt_pos_etd,strlen(ref),CNV_1_str);
                  //Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_etd,strlen(ref),CNV_2_str);
                }
                
              }

            }
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }
        }
        else{
          while (vcfpos < start_pos){
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }
        }
        start_pos += dist + 1;
      }
    }
    chr_no++;
  }
  //bcf_hdr_destroy(test_header);
  //bcf_destroy(test_record); 
  //bcf_close(test_vcf);  
  return; 
}


void VCF3(const char* fastafile,double cov=1.0){  
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
  test_vcf = vcf_open("/home/wql443/WP1/data/chr14_full.vcf", "r");
  test_header = bcf_hdr_read(test_vcf);

  while (chr_no < faidx_nseq(ref)){
    // The fasta file from the references
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);

    int start_pos =1;// 19108144; 
    int end_pos = name_len;19111144;  

    //14	19110144	.	GTGGCTGCAGCCA

    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    //std::cout << "ts" << chrno << std::endl;
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos = test_record->pos+1;

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
        else if (start_pos <= vcfpos){
          std::cout << "------------ NEW REGION -------" << std::endl;
          char nt[] = "tT";
          
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int length = strlen(sequence);
          
          int alt_bool = 0;

          //extend the sequence
          int alt_pos_ext=0;
          int alt2_pos_ext=0;
          
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
              if (vcfpos-start_pos > 0){alt_pos = vcfpos-start_pos-1;}
              else if(vcfpos-start_pos == 0){alt_pos = vcfpos-start_pos;} //to ensure variation in the first basepair..

              ref = test_record->d.allele[0];
              alt = test_record->d.allele[1];

              alt_no = test_record->d.info->len;

              std::string refstr(ref);
              std::string altstr(alt);

              int refsize = (int) refstr.size();
              int altsize = (int) altstr.size();

              if (alt_no < 2){
                //single alternative
                ref_init_len = refsize;
                std::cout << " single " << refstr << " alt " << altstr << "size "<< altsize << std::endl;
                std::cout << "alt " << alt_pos << " ext " << alt_pos_ext << " new pos " << alt_pos + alt_pos_ext <<" len " << Seq_str.length() << std::endl;
                
                // single GTGTGTGGA && (alt_pos+alt_pos_ext)>ref_init_len
                if ((alt_pos+alt_pos_ext)<0 && (alt_pos+alt_pos_ext)>ref_init_len){ //&& (alt_pos+alt_pos_ext) > start_pos
                  std::cout << "if " << std::endl;
                  //vcfstruct = bcf_read(test_vcf, test_header, test_record);
                  //vcfpos = test_record->pos+1;
                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                }
                else if ((alt_pos+alt_pos_ext)>=0){
                  std::cout << "eles if" << std::endl;
                  std::cout << Seq_str << Seq_str.length() << std::endl;
                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                  std::cout << Seq_str << Seq_str.length() << std::endl;
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                }
                else
                {
                  std::cout << "else "<< std::endl;
                  Seq_str = Seq_str;
                  Seq_str_alt = Seq_str_alt;
                  std::cout << Seq_str << Seq_str.length() << std::endl;
                }
                

                //alt_pos_ext += 1;
                if(refsize == altsize){
                  //substituion
                  std::cout << "sub" << std::endl;
                  alt_pos_ext += (refsize-altsize);
                }
                else if (refsize > altsize){
                  //deletion
                  std::cout << "del" << std::endl;
                  alt_pos_ext -= (refsize-altsize);
                }
                else if (altsize > refsize){
                  //insertions
                  alt_pos_ext += std::abs(refsize-altsize);
                }
                
                std::cout << "alt pos test "<< alt_pos_ext << std::endl;
                /*if(refsize > altsize){
                // deletions of ref to alternative
                std::cout << "if redsize" << std::endl;
                std::cout << "alternative " << alt_pos << std::endl;
                std::cout << alt_pos_ext << " ref " << refsize << " alt " << altsize << std::endl;
                alt_pos_ext = alt_pos_ext - (refsize-altsize);
                std::cout << alt_pos_ext << " ref2 " << refsize << " alt2 " << altsize << std::endl;
                //&& (alt_pos + alt_pos_ext) < 0
                if (alt_pos_ext < 0){
                  std::cout << "lol" << std::endl;
                  alt_pos_ext = 0;
                }
                else
                {
                  alt_pos_ext = alt_pos_ext - (refsize-altsize);
                }*/

              }
              else{
                //multiple alternative
                alt_bool += 1;

                alt2 = test_record->d.allele[alt_no];
                std::string altstr2(alt2);
                int alt2size = (int) altstr2.size();  
                std::cout <<"-------------------" << std::endl;
                std::cout << "multiple ref" << ref << refsize << std::endl;
                std::cout << "multiple " << alt << " lol "<< alt2 << " size " << alt2size << std::endl;

                if (alt[0] == '<'){
                  std::cout << "alternative if"<<std::endl;
                  //14	19112155	DUP_gs_CNV_14_19112155_19126473	T	<CN2>,<CN3>
                  //removes <CN from all the copy number alternatives and extracts the number //it might give problems with X < 9 CNV10 etc..
                  altstr.erase(0,3).pop_back();  
                  altstr2.erase(0,3).pop_back();
                  int CNV_1 = stoi(altstr);
                  int CNV_2 = stoi(altstr2);
                  
                  //creates the copy number string
                  std::string CNV_1_str(CNV_1, ref[0]); 
                  std::string CNV_2_str(CNV_2, ref[0]); 

                  if(refsize > CNV_2){alt2_pos_ext = alt2_pos_ext - std::abs(refsize-CNV_2);}
                  else if(refsize <= CNV_2){alt2_pos_ext = alt2_pos_ext + std::abs(refsize-CNV_2);}

                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),CNV_1_str);

                  //For some reason it fucks up if i say alt_pos+alt2_pos_ext 
                  Seq_str_alt = Seq_str_alt.replace(alt_pos,strlen(ref),CNV_2_str);
                }
                else
                {  
                  //multiple different unique alternatives
                  std::cout << "alternative else "<< std::endl;
                  std::cout << "alter" << alt_pos_ext << " - " << alt2_pos_ext << std::endl;
                  if(refsize > alt2size){alt2_pos_ext = alt2_pos_ext - std::abs(refsize-alt2size);}
                  else if(refsize <= alt2size){alt2_pos_ext = alt2_pos_ext + std::abs(refsize-alt2size);}
                  std::cout << "alter2" << alt_pos_ext << " - " << alt2_pos_ext << std::endl;        
                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                  std::cout << "first succes "<< std::endl;
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt2_pos_ext,strlen(ref),altstr2);
                }
              }
            
              // extending the positions to ensure correct index for the sequences
              /*if(refsize > altsize){
                // deletions of ref to alternative
                std::cout << "if redsize" << std::endl;
                std::cout << "alternative " << alt_pos << std::endl;
                std::cout << alt_pos_ext << " ref " << refsize << " alt " << altsize << std::endl;
                alt_pos_ext = alt_pos_ext - (refsize-altsize);
                std::cout << alt_pos_ext << " ref2 " << refsize << " alt2 " << altsize << std::endl;
                //&& (alt_pos + alt_pos_ext) < 0
                if (alt_pos_ext < 0){
                  std::cout << "lol" << std::endl;
                  alt_pos_ext = 0;
                }
                else
                {
                  alt_pos_ext = alt_pos_ext - (refsize-altsize);
                }*/
                
                /*if ((alt_pos_ext - std::abs(refsize-altsize))<0)
                {//Ensure you know what youre doing here..
                  alt_pos_ext = 0;
                }
                else
                {
                  alt_pos_ext = alt_pos_ext - std::abs(refsize-altsize);
                }
                std::cout << alt_pos_ext << " changed " << std::abs(refsize-altsize) << std::endl;
                }

              else if(refsize <= altsize){alt_pos_ext += std::abs(refsize-altsize);}*/
            }
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }

          if (alt_bool == 0){
            //std::cout << "vcf_seq" << std::endl << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
            outfa << Seq_str << std::endl;
          }
          else{
            //for 2 alternatives create two seperate sequences, NB! improve for higher number of alternatives.
            //std::cout << "vcf_seq" << std::endl << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
            outfa << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_alt2_len:" << length << "_" << Seq_str_alt.length() << std::endl;
            outfa << Seq_str_alt << std::endl;
          }
        }

        else{
          while (vcfpos < start_pos){
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }
        }
        start_pos += dist + 1;
      }
    }
    chr_no++;
  }
  //bcf_hdr_destroy(test_header);
  //bcf_destroy(test_record); 
  //bcf_close(test_vcf);  
  return; 
}

void VCF4(const char* fastafile,double cov=1.0){  
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

    int start_pos =20250390;//20250390;// 19108144; 
    int end_pos = 20255390;//20255390; //19111144;  

    //14	19110144	.	GTGGCTGCAGCCA
    //14	19291899	.	C	CTGTT	100 //19291650; //19291900; -> CTGTTG (WORKS)
    //start 20265509 vcf 20265521 end 20265563
    //ref  CATTTGCCAGGTA alt C
    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    //std::cout << "ts" << chrno << std::endl;
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos = test_record->pos+1;

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
        else if (start_pos <= vcfpos){
          std::cout << "------------ NEW REGION -------" << std::endl;
          char nt[] = "tT";
          
          //Deamin_char(sequence,nt,rand_len);
          std::string Seq_str(sequence);
          std::string Seq_str_alt(sequence);
          int length = strlen(sequence);
          
          int alt_bool = 0;

          //extend the sequence
          int alt_pos_ext=0;
          int alt2_pos_ext=0;
          
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
              if (vcfpos-start_pos > 0){alt_pos = vcfpos-start_pos-1;}
              else if(vcfpos-start_pos == 0){alt_pos = vcfpos-start_pos;} //to ensure variation in the first basepair..

              ref = test_record->d.allele[0];
              alt = test_record->d.allele[1];
              std::cout << "ref  " << ref  << " alt "<<  alt << std::endl;

              alt_no = test_record->d.info->len;

              std::string refstr(ref);
              std::string altstr(alt);

              int refsize = (int) refstr.size();
              int altsize = (int) altstr.size();
              int new_alt_pos = alt_pos+alt_pos_ext;
              if (alt_no < 2){
                //single alternative
                ref_init_len = refsize;
                std::cout << " single " << refstr << " alt " << altstr << "size "<< altsize << std::endl;
                std::cout << "alt " << alt_pos << " ext " << alt_pos_ext << " new pos " << new_alt_pos <<" len " << Seq_str.length() << std::endl;
                
                // single GTGTGTGGA && (alt_pos+alt_pos_ext)>ref_init_len
                if ((new_alt_pos)<0 && (new_alt_pos)>ref_init_len){ //&& (alt_pos+alt_pos_ext) > start_pos
                  std::cout << "if " << std::endl;
                  //vcfstruct = bcf_read(test_vcf, test_header, test_record);
                  //vcfpos = test_record->pos+1;
                  Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                  Seq_str_alt = Seq_str_alt.replace(new_alt_pos,strlen(ref),altstr);
                }
                else if (0<=new_alt_pos && new_alt_pos<Seq_str_alt.length()){
                  std::cout << "eles if" << std::endl;
                  std::cout << "ref before " << Seq_str << Seq_str.length() << std::endl;
                  Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                  std::cout << "ref  after " << Seq_str << Seq_str.length() << std::endl;
                  std::cout << "alt before " << Seq_str_alt << Seq_str_alt.length() << std::endl;
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                  std::cout << "ref  after " << Seq_str_alt << Seq_str_alt.length() << std::endl;

                }
                else
                {
                  //std::cout << "else "<< std::endl;
                  Seq_str = Seq_str;
                  Seq_str_alt = Seq_str_alt;
                  //std::cout << Seq_str << Seq_str.length() << std::endl;
                }
                

                //alt_pos_ext += 1;
                if(refsize == altsize){
                  //substituion
                  //std::cout << "sub" << std::endl;
                  alt_pos_ext += (refsize-altsize);
                }
                else if (refsize > altsize){
                  //deletion
                  //std::cout << "del" << std::endl;
                  alt_pos_ext -= (refsize-altsize);
                }
                else if (altsize > refsize){
                  //insertions
                  alt_pos_ext += std::abs(refsize-altsize);
                }
              }
              else{
                //multiple alternative
                alt_bool += 1;

                alt2 = test_record->d.allele[alt_no];
                std::string altstr2(alt2);
                int alt2size = (int) altstr2.size();  
                std::cout << "mulitple alternatives "<< std::endl;
                std::cout << "seq length" << Seq_str.length() << std::endl;
                std::cout << "ref " << ref << " size " << refsize << std::endl;
                std::cout << "alt1 " << alt << " alt2 "<< alt2 << std::endl;
                
                if (alt[0] == '<'){
                  //std::cout << "alternative if"<<std::endl;
                  //14	19112155	DUP_gs_CNV_14_19112155_19126473	T	<CN2>,<CN3>
                  //removes <CN from all the copy number alternatives and extracts the number //it might give problems with X < 9 CNV10 etc..
                  altstr.erase(0,3).pop_back();  
                  altstr2.erase(0,3).pop_back();
                  int CNV_1 = stoi(altstr);
                  int CNV_2 = stoi(altstr2);
                  
                  //creates the copy number string
                  std::string CNV_1_str(CNV_1, ref[0]); 
                  std::string CNV_2_str(CNV_2, ref[0]); 

                  if(refsize > CNV_2){alt2_pos_ext = alt2_pos_ext - std::abs(refsize-CNV_2);}
                  else if(refsize <= CNV_2){alt2_pos_ext = alt2_pos_ext + std::abs(refsize-CNV_2);}

                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),CNV_1_str);

                  //For some reason it fucks up if i say alt_pos+alt2_pos_ext 
                  Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_ext,strlen(ref),CNV_2_str);
                }
                else
                {  
                  // if the alternative alleles are on the first basepair the index are affected as its normally both insertion and deltion of the ref

                  /*if ((alt2_pos_ext+vcfpos) < start_pos || (alt2_pos_ext+vcfpos) < start_pos)
                  {

                  }*/
                  
                  if (vcfpos-start_pos == 0) //keep index in seq for variation if first basepar
                  {     
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr);
                    Seq_str_alt = Seq_str_alt.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else if ((vcfpos-refsize)-start_pos < 0) //variation in the first basepairs with one alternative -> deletions removing index 0
                  {
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr);
                    Seq_str_alt = Seq_str_alt.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else
                  {
                    //multiple different unique alternatives
                    std::cout << "alternative else "<< std::endl;
                    std::cout << "alt pos "<< alt_pos << std::endl;
                    std::cout << Seq_str << std::endl;
                    std::cout << "alter:" << alt_pos_ext << " alter 2 " << alt2_pos_ext << std::endl;
                    if(refsize > alt2size){alt2_pos_ext = alt_pos_ext - std::abs(refsize-alt2size);}
                    else if(refsize <= alt2size){alt2_pos_ext = alt2_pos_ext + std::abs(refsize-alt2size);}
                  
                    std::cout << "alter2:" << alt_pos_ext << " alter 2 " << alt2_pos_ext << std::endl;        
                    Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr);
                    std::cout << "first succes "<< std::endl;
                    std::cout << Seq_str << std::endl;
                    Seq_str_alt = Seq_str_alt.replace(alt_pos+alt_pos_ext,strlen(ref),altstr2); //alt2_pos_ext
                    std::cout << "second succes "<< std::endl;
                    std::cout << Seq_str_alt << std::endl;
                  }
                  //14	22015934	rs73581477	A	G	100
                }
              }
            std::cout << " --- end of loop --- " << std::endl;
            }
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }

          if (alt_bool == 0){
            //std::cout << "vcf_seq" << std::endl << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
            outfa << Seq_str << std::endl;
          }
          else{
            //for 2 alternatives create two seperate sequences, NB! improve for higher number of alternatives.
            //std::cout << "vcf_seq" << std::endl << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
            outfa << Seq_str << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_alt2_len:" << length << "_" << Seq_str_alt.length() << std::endl;
            outfa << Seq_str_alt << std::endl;
          }
        }

        else{
          while (vcfpos < start_pos){
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }
        }
        start_pos += dist + 1;
      }
    }
    chr_no++;
  }
  //bcf_hdr_destroy(test_header);
  //bcf_destroy(test_record); 
  //bcf_close(test_vcf);  
  return; 
}


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

    int start_pos =1; //20250390;//20250390;// 19108144; 
    int end_pos = name_len;//;20255390;//20255390; //19111144;  

    //14	19110144	.	GTGGCTGCAGCCA
    //14	19291899	.	C	CTGTT	100 //19291650; //19291900; -> CTGTTG (WORKS)
    //start 20265509 vcf 20265521 end 20265563
    //ref  CATTTGCCAGGTA alt C
    //create the proper chromosome name for VCF format
    int vcfstruct = bcf_read(test_vcf, test_header, test_record);
    char vcfchr[4096];   // array to hold the result.
    const char *chr = "chr";
    const char *chrno = bcf_hdr_id2name(test_header, test_record->rid);
    strcpy(vcfchr,chr);
    strcat(vcfchr,chrno);

    //ensure same chromosome from fasta and VCF
    if (std::strcmp(vcfchr,name) == 0){
      int vcfpos = test_record->pos+1;

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
        else if (start_pos <= vcfpos){
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
              if (vcfpos-start_pos > 0){alt_pos = vcfpos-start_pos-1;}
              else if(vcfpos-start_pos == 0){alt_pos = vcfpos-start_pos;} //to ensure variation in the first basepair..

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
                
                //(new_alt_pos)<0 updated position for deletion in beginning
                //(new_alt_pos)>ref_init_len  updated position is higher than the length of the original ref allele?
                if ((new_alt_pos)<0 && (new_alt_pos)>start_pos){
                  std::cout << "only if" << std::endl;
                  Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                }
                //variation within the "middle"
                else if (0<=new_alt_pos && new_alt_pos<Seq_str_alt.length()){
                  std::cout << "else if "<< std::endl;
                  std::cout << "alt pos" << new_alt_pos << std::endl;
                  Seq_str = Seq_str.replace(new_alt_pos,strlen(ref),altstr);
                }
                //No variation
                else
                {
                  std::cout << "ELSE" << std::endl;
                  std::cout << "alt pos" << new_alt_pos << std::endl;
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
                //selecting one random of the multiple random alleles
                int rand_alt = rand() % alt_no + 1;
                std::cout << "RANDOM" << rand_alt << " akt " << alt_no << std::endl;

                alt2 = test_record->d.allele[rand_alt];
                std::string altstr2(alt2);
                int alt2size = (int) altstr2.size();  
                std::cout << "mulitple alternatives "<< std::endl;
                std::cout << "seq length" << Seq_str.length() << std::endl;
                std::cout << "ref " << ref << " size " << refsize << std::endl;
                std::cout << "alt1 " << alt << " alt2 "<< alt2 << std::endl;
                
                if (alt[0] == '<'){
                  //std::cout << "alternative if"<<std::endl;
                  //14	19112155	DUP_gs_CNV_14_19112155_19126473	T	<CN2>,<CN3>
                  //removes <CN from all the copy number alternatives and extracts the number //it might give problems with X < 9 CNV10 etc..
                  altstr.erase(0,3).pop_back();  
                  altstr2.erase(0,3).pop_back();
                  int CNV_1 = stoi(altstr);
                  int CNV_2 = stoi(altstr2);
                  
                  //creates the copy number string
                  std::string CNV_1_str(CNV_1, ref[0]); 
                  std::string CNV_2_str(CNV_2, ref[0]); 

                  Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),CNV_1_str);

                }
                else
                {  
                  // if the alternative alleles are on the first basepair the index are affected as its normally both insertion and deltion of the ref
                  if (vcfpos-start_pos == 0) //keep index in seq for variation if first basepar
                  { 
                    std::cout <<  "if "<< std::endl;    
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else if ((vcfpos-refsize)-start_pos < 0) //variation in the first basepairs with one alternative -> deletions removing index 0
                  {
                    std::cout <<  "else if "<< std::endl;
                    Seq_str = Seq_str.replace(alt_pos,strlen(ref),altstr2);
                  }
                  else
                  {
                    std::cout <<  "else "<< std::endl;
                    //multiple different unique alternatives
                    Seq_str = Seq_str.replace(alt_pos+alt_pos_ext,strlen(ref),altstr2);
                  }
                  //14	22015934	rs73581477	A	G	100
                  if(refsize == alt2size){alt_pos_ext += 0;}
                  else if (refsize > alt2size){alt_pos_ext -= (refsize-alt2size);}
                  else if (alt2size > refsize){alt_pos_ext += (alt2size-refsize);}
                }
              }
            std::cout << " --- end of loop --- " << std::endl;
            }
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
          }

          outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_len:" << length << "_" << Seq_str.length() << std::endl;
          outfa << Seq_str << std::endl;
          
        }

        else{
          while (vcfpos < start_pos){
            //iterating through the vcf file
            vcfstruct = bcf_read(test_vcf, test_header, test_record);
            vcfpos = test_record->pos+1;
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
