#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <math.h>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"
#include "mrand.h"

#define LENS 4096
#define MAXBINS 100
unsigned char nuc2int[255];

pthread_mutex_t Fq_write_mutex;

int MacroRandType = 0;

struct Parsarg_for_Sampling_thread{
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  char *genome;
  int *chr_idx_array;
  int chr_no;
  int threadno;
  size_t *size_cumm;
  const char **names;

  int* FragLen;
  double* FragFreq;
  int No_Len_Val;

  char *NtQual_r1;
  char *NtQual_r2;
  ransampl_ws ***QualDist_r1;
  ransampl_ws ***QualDist_r2;

  double* MisMatch;
  const char* SubFlag;
  int MisLength;

  double *NtErr_r1;
  double *NtErr_r2;

  int threadseed;
  size_t reads;
  
  BGZF *bgzf_fp1;
  BGZF *bgzf_fp2;

  samFile *SAMout;
  sam_hdr_t *SAMHeader;
  bam1_t **list_of_reads;
  int l;
  int m;

  const char* Adapter_flag;
  const char* Adapter_1;
  const char* Adapter_2;

  const char* Briggs_flag;
  float *BriggsParam;
  
  const char* SizeFile;
  const char* SizeFileFlag;
  int FixedSize;

  int readcycle;
  
  const char* OutputFormat;
  const char* SeqType;
  const char* QualFlag;

  char ErrorFlag;
  char PolyNt;

  int RandMacro;
};
      
void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  int MacroRandType;
  MacroRandType = struct_obj -> RandMacro; //LLLLLL
  //fprintf(stderr,"RANDOM VALUE v1 %d \n",MacroRandType);
  //fprintf(stderr,"RANDOM VALUE2 %d \n",MacroRandType);
  // creating random objects for all distributions.
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  size_t genome_len = strlen(struct_obj->genome);

  //int MacroRandType = struct_obj -> RandMacro; //LLLLLL

  mrand_t *drand_alloc = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_nt = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_nt_adapt = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_briggs = mrand_alloc(MacroRandType,loc_seed);
  //fprintf(stderr,"Macro type %d \t Seed val %d\n",MacroRandType,loc_seed);
   
  // sequence reads, original, modified, with adapters, pe
  char seq_r1[1024] = {0};
  char seq_r1_mod[1024] = {0};

  char read[1024] = {0};
  char readadapt[1024] = {0};

  char seq_r2[1024] = {0};
  char seq_r2_mod[1024] = {0};
  char read2[1024] = {0};
  char readadapt2[1024] = {0};
  
  char read_rc_sam1[1024] = {0};
  char readadapt_rc_sam1[1024] = {0};
  char read_rc_sam2[1024] = {0};
  char readadapt_rc_sam2[1024] = {0};
  char readadapt_err1[1024] = {0};
  char readadapt_err2[1024] = {0};

  int reads = struct_obj -> reads;

  size_t rand_start;

  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};

  int localread = 0;
  int iter = 0;
  int current_reads_atom = 0;
  int readsizelimit;
   
  while (current_reads_atom < reads){
    /*double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX); //random between 0 and 1
    double rand_val2 = rand_val*RAND_MAX; // between 0 and maximum 
    double rand_val3 = myrand((unsigned int) rand_val2); //((double) rand_r(&test)/ RAND_MAX);// 
    rand_start = rand_val3 * (genome_len-300); //genome_len-100000;*/
    double rand_val = mrand_pop(drand_alloc);
    rand_start = rand_val * (genome_len-300)+1; //genome_len-100000;
    //fprintf(stderr,"random start2 %zu\n-----------\n",rand_start2);
    // Fragment length creation
    int fraglength;
    if (struct_obj->No_Len_Val != -1){
      // random start and length are not dependent on the same rand val
      //fprintf(stderr,"Legth%f \n",mrand_pop(drand_alloc));
      int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val);
      fraglength =  struct_obj->FragLen[lengthbin];
      //fprintf(stderr,"rand val %f \t fraglent %d\n",rand_val,fraglength);
    }
    else{
      //fprintf(stderr,"FIXED LENGTH \n");
      fraglength = struct_obj->FixedSize;
    } 
    if(strcasecmp("false",struct_obj->QualFlag)==0){
      readsizelimit = fraglength;
    }
    else{
      readsizelimit = struct_obj->readcycle;
    }
    //we need a new random value here otherwise the random ID would be the same for 
    //based on the fragment lengths
    double rand_val_id = mrand_pop(drand_alloc);
    //int rand_id = (rand_val * fraglength-1); //100
    int rand_id = (rand_val_id * fraglength-1); //100
    //fprintf(stderr,"Random val %f \t id %d\n",rand_val_id,rand_id);
    
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}
    //fprintf(stderr,"chromosome index 2 %d\n",chr_idx);
    if (fraglength > readsizelimit){strncpy(seq_r1,struct_obj->genome+rand_start-1,readsizelimit);}   // case 1
    else {strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    
    if(strcasecmp("PE",struct_obj->SeqType)==0){
      if (fraglength > readsizelimit){

        strncpy(seq_r2,struct_obj->genome+rand_start+fraglength-1-readsizelimit,readsizelimit);
        } // case 1
      else {strncpy(seq_r2,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    }

    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N');
    pch2 = strrchr(seq_r1,'N');

    int seqlen = strlen(seq_r1);
    //std::cout << seq_r1 << std::endl;
    //fprintf(stderr,"random_start %zu and chromosome index %d with chromosome size %zu and new size %zu\n",rand_start,chr_idx,struct_obj->size_cumm[chr_idx],struct_obj->size_cumm[chr_idx+1]);
    //if((rand_start-struct_obj->size_cumm[chr_idx])<seqlen){fprintf(stderr,"LOOOOOOORT %zu \t %zu \t %d \n",rand_start,struct_obj->size_cumm[1],seqlen); exit(0);}

    size_t n_cigar;const uint32_t *cigar;const uint32_t *cigar2; uint16_t flag; uint16_t flag2;
    uint32_t cigar_bitstring = bam_cigar_gen(seqlen, BAM_CMATCH);

    if ((int )(pch-seq_r1+1) == 1 || (int)(pch2-seq_r1+1)  == seqlen){
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);}
    else if ((rand_start-struct_obj->size_cumm[chr_idx+1])<fraglength){ //(rand_start-struct_obj->size_cumm[chr_idx])<fraglength
      // perhaps readsizelimit or seqlen
      /*fprintf(stderr,"------------------------\n");
      std::cout << "RANDOM START LOOP" << std::endl;
      std::cout << seq_r1 << std::endl;
      fprintf(stderr,"random_start %zu and chromosome index %d with chromosome size %zu and new size %zu\n",rand_start,chr_idx,struct_obj->size_cumm[chr_idx],struct_obj->size_cumm[chr_idx+1]);*/
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);
    }
    else{
      // then all the same start pos would have the same strand no matter the chromosome?
      int strand = (int) rand_start%2;//(int) rand_r(&loc_seed)%2;//1;//rand() % 2;
      //fprintf(stderr,"STRAND %d\n",strand);
      //fprintf(stderr,"STRAND %d\t%d\t%f\t%d\n",rand_r(&loc_seed),(int) rand_r(&loc_seed)%2,mrand_pop(drand_alloc),rand_start%2);
      if (struct_obj->SAMout){
        // in bam files the reads are all aligning to the forward strand, but the flags identify read property
        if (strcasecmp("SE",struct_obj->SeqType)==0){
          if (strand == 0){flag = 0;} // se read forward strand
          else{flag = 16;
            //fprintf(stderr,"--------------\nFLAG IS %d\n ORIG\n%s \n",flag,seq_r1);
            // Simulate reads from the negative strand but change the orientation to 5' to 3' fragment
            DNA_complement(seq_r1);
            reverseChar(seq_r1,strlen(seq_r1));
            //fprintf(stderr,"FIRST REV COMP\n%s \n",seq_r1);
          }
        }
        if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 0){
          flag = 97;  //Paired, mate reverse, first
          flag2 = 145; // Paired, reverse strand, second
          DNA_complement(seq_r2);
          reverseChar(seq_r2,strlen(seq_r2));
        }
        else if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 1){
          flag = 81; // Paired, reverse strand, first
          flag2 = 161; // Paired, mate reverse, second
          DNA_complement(seq_r1);
          reverseChar(seq_r1,strlen(seq_r1));
        }

        // To add deamination we need the fragments to emulate both forward and reverse strand with the real orientation (hence the reverse complement earlier)
        if(strcasecmp(struct_obj->Briggs_flag,"true")==0){
          SimBriggsModel(seq_r1, seq_r1_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
          
          strncpy(seq_r1, seq_r1_mod, sizeof(seq_r1));
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            SimBriggsModel(seq_r2, seq_r2_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
            strncpy(seq_r2, seq_r2_mod, sizeof(seq_r2));
          }
        }

        // Similar if we use deamination file
        if(strcasecmp("true",struct_obj->SubFlag)==0){
          Deam_File(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            Deam_File(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          }
        }

        //I chose not to add sequencing errors here, since we need the adapter to be added first.

        // HERE WE USE SEQUENCING ERRORS TO CREATE SUBSTITUTIONS.
        // ADD SEQUENCING ERROR HERE IN ORDER TO GET CORRECT ORIENATION OF THE ADDED NUCLEOTIDE SUBSTITUTIONS

        // Due to real life empirical sam files and its specifications we need all single-end data to be from the positive strand, as such we need to 
        // change the orientation and sequence for the reads from the reverse strand back to the forward strand.
        
        /*if (flag == 16 || flag == 81){
          fprintf(stderr,"INSIDE DNA COMPLEMENT LOOP\n------------------\n");
          DNA_complement(seq_r1);reverseChar(seq_r1);}
        else if (flag == 97){DNA_complement(seq_r2);reverseChar(seq_r2);}*/  
      }
      else{
        // in fasta and fastq the sequences need to be on forward or reverse strand, i.e we need reverse complementary        
        if (strcasecmp("SE",struct_obj->SeqType)==0 && strand == 1){
          DNA_complement(seq_r1);
          reverseChar(seq_r1,strlen(seq_r1));
        }

        if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 0){
          DNA_complement(seq_r2);
          reverseChar(seq_r2,strlen(seq_r2));
        }
        else if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 1){
          DNA_complement(seq_r1);
          reverseChar(seq_r1,strlen(seq_r1));
        }

        // Adding PMD after strand is selected, as to not influence the symmetry of the PMD
        if(strcasecmp(struct_obj->Briggs_flag,"true")==0){
          SimBriggsModel(seq_r1, seq_r1_mod, strlen(seq_r1),struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
          strncpy(seq_r1, seq_r1_mod, sizeof(seq_r1));

          if (strcasecmp("PE",struct_obj->SeqType)==0){
            SimBriggsModel(seq_r2, seq_r2_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
            strncpy(seq_r2, seq_r2_mod, sizeof(seq_r2));
          }
        }

        if(strcasecmp("true",struct_obj->SubFlag)==0){
          //fprintf(stderr,"mismatch matrix %f \t %f \t %f \t %f \t %f \t %f \n",struct_obj->MisMatch[0],struct_obj->MisMatch[10],struct_obj->MisMatch[14],struct_obj->MisMatch[15],struct_obj->MisMatch[30],struct_obj->MisMatch[45]);
          //fprintf(stderr,"SEQUENCE %s\n",seq_r1);
          Deam_File(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          // IS THIS THE WAY TO DO IT?
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            Deam_File(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          }
          //fprintf(stderr,"SEQUENCE %s\n",seq_r1);
        }
      }
      
      char READ_ID[1024]; int read_id_length;
      read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%zu-%zu_length:%d", struct_obj->threadno, rand_id,strand,struct_obj->names[chr_idx],
        rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],fraglength);
      
      if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
        strcpy(read, seq_r1);
        strcat(read,struct_obj->Adapter_1);

        if(strcasecmp("false",struct_obj->QualFlag)==0){
          strncpy(readadapt, read, strlen(read));
          readsizelimit = strlen(read);
        }
        else{strncpy(readadapt, read, readsizelimit);}

        if (strcasecmp("PE",struct_obj->SeqType)==0){
          strcpy(read2, seq_r2);
          strcat(read2,struct_obj->Adapter_2);
          if(strcasecmp("false",struct_obj->QualFlag)==0){
            strncpy(readadapt2, read2, strlen(read2));
            readsizelimit = strlen(read2);
          }
          else{strncpy(readadapt2, read2, readsizelimit);}
        }

        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,readadapt2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          for(long unsigned int p = 0;p<strlen(readadapt);p++){
            //double rand_val = mrand_pop(drand_alloc);
            double dtemp1;double dtemp2;
            dtemp1 = mrand_pop(drand_alloc_nt_adapt);
            dtemp2 = mrand_pop(drand_alloc_nt_adapt);
            //fprintf(stderr,"temp lol 1 %f \t temp 2 %f\n",dtemp1,dtemp2);

            int base = readadapt[p];
            int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
            qual_r1[p] = struct_obj->NtQual_r1[qscore];
            
            if (struct_obj->ErrorFlag == 'T'){
              double dtemp3;double dtemp4;
              dtemp3 = mrand_pop(drand_alloc_nt_adapt);
              dtemp4 = mrand_pop(drand_alloc_nt_adapt);
              if (dtemp3 < struct_obj->NtErr_r1[qscore]){ErrorSub(dtemp4,readadapt,p);}
            }
          }

          if (struct_obj->PolyNt != 'F'){
            char PolyChar[1024] = "\0"; char QualPoly[1024]= "\0";
            memset(PolyChar,struct_obj->PolyNt, readsizelimit);
            memset(QualPoly,'!', readsizelimit);
            strncpy(PolyChar, readadapt, strlen(readadapt));
            strncpy(QualPoly, qual_r1, strlen(qual_r1));
            ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,PolyChar,QualPoly);
          }
          else{ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,readadapt,qual_r1);}
          
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            for(long unsigned int p = 0;p<strlen(readadapt2);p++){
              double dtemp1;double dtemp2;
              dtemp1 = mrand_pop(drand_alloc_nt_adapt);
              dtemp2 = mrand_pop(drand_alloc_nt_adapt);

              int base = readadapt2[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
              qual_r2[p] = struct_obj->NtQual_r2[qscore];
              if (struct_obj->ErrorFlag == 'T'){
                double dtemp3;double dtemp4;
                dtemp3 = mrand_pop(drand_alloc_nt_adapt);
                dtemp4 = mrand_pop(drand_alloc_nt_adapt);
                if (dtemp3 < struct_obj->NtErr_r2[qscore]){ErrorSub(dtemp4,readadapt2,p);}
              }
            }
            
            if (struct_obj->PolyNt != 'F'){
              char PolyChar[1024]; char QualPoly[1024];
              memset(PolyChar,struct_obj->PolyNt, readsizelimit);
              memset(QualPoly,'!', readsizelimit);
              strncpy(PolyChar, readadapt2, strlen(readadapt2));
              strncpy(QualPoly, qual_r2, strlen(qual_r2));
              ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,PolyChar,QualPoly);
            }
            else{ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,readadapt2,qual_r2);}
          }
        }
        if (struct_obj->SAMout){
          if(strcasecmp("true",struct_obj->QualFlag)==0){
            // if we actually want the nucleotide quality string in the sam output
            for(long unsigned int p = 0;p<strlen(readadapt);p++){
              double dtemp1;double dtemp2;
              dtemp1 = mrand_pop(drand_alloc_nt_adapt);
              dtemp2 = mrand_pop(drand_alloc_nt_adapt);

              int base = readadapt[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
              qual_r1[p] = struct_obj->NtQual_r1[qscore];
              if (struct_obj->ErrorFlag == 'T'){
                double dtemp3;double dtemp4;
                dtemp3 = mrand_pop(drand_alloc_nt_adapt);
                dtemp4 = mrand_pop(drand_alloc_nt_adapt);
                if (dtemp3 < struct_obj->NtErr_r1[qscore]){ErrorSub(dtemp4,readadapt,p);}
              }
            }
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              for(long unsigned int p = 0;p<strlen(readadapt2);p++){
                double dtemp1;double dtemp2;
                dtemp1 = mrand_pop(drand_alloc_nt_adapt);
                dtemp2 = mrand_pop(drand_alloc_nt_adapt);

                int base = readadapt2[p];
                int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
                qual_r2[p] = struct_obj->NtQual_r2[qscore];
                
                if (struct_obj->ErrorFlag == 'T'){
                  double dtemp3;double dtemp4;
                  dtemp3 = mrand_pop(drand_alloc_nt_adapt);
                  dtemp4 = mrand_pop(drand_alloc_nt_adapt);

                  if (dtemp3 < struct_obj->NtErr_r2[qscore]){ErrorSub(dtemp4,readadapt2,p);}
                }
              }
            }
          }

          if (flag == 0){sprintf(readadapt_rc_sam1, "%s",readadapt);}
          else if (flag == 16 || flag == 81){
            
            // 16 se reverse strand // 81 first in pair reverse strand
            // we cannot simply copy the adapter sequence, since sequencing error is possible, as such
            // readadapt_err is the extracted sequence after Suberr 
            sprintf(readadapt_err1, "%*s", strlen(readadapt)-strlen(seq_r1), readadapt+strlen(seq_r1));
            sprintf(read_rc_sam1, "%.*s", strlen(seq_r1), readadapt);
            DNA_complement(read_rc_sam1);reverseChar(read_rc_sam1,strlen(read_rc_sam1));
            sprintf(readadapt_rc_sam1, "%s%s", read_rc_sam1,readadapt_err1);
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              // the mate to flag 81 will be the second in pair on the forward strand - just to esnrue identical names
              sprintf(readadapt_rc_sam2, "%s",readadapt2);
            } 
          }
          else if (flag == 97){
            // 97 first in pair forward strand only for Paired end.
            sprintf(readadapt_rc_sam1, "%s",readadapt); //to ensure all the reads have identifcal names for the sam output

            // the mate is reverse
            sprintf(readadapt_err2, "%*s", strlen(readadapt2)-strlen(seq_r2), readadapt2+strlen(seq_r2));
            sprintf(read_rc_sam2, "%.*s", strlen(seq_r2), readadapt2); // copy the sequence 2 from the sequence 2 + adapter
            DNA_complement(read_rc_sam2);reverseChar(read_rc_sam2,strlen(read_rc_sam1)); // reverse complement sequence 2 
            sprintf(readadapt_rc_sam2, "%s%s", read_rc_sam2,readadapt_err2); // create rev comp sequence 2 + adapter
          }

          if (struct_obj->PolyNt != 'F'){
            char PolyChar[1024] = "\0";
            memset(PolyChar,struct_obj->PolyNt, readsizelimit);
            strncpy(PolyChar, readadapt_rc_sam1, strlen(readadapt_rc_sam1));
            ksprintf(struct_obj->fqresult_r1,"%s",PolyChar);
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              char PolyChar[1024] = "\0";
              memset(PolyChar,struct_obj->PolyNt, readsizelimit);
              strncpy(PolyChar, readadapt_rc_sam2, strlen(readadapt_rc_sam2));
              ksprintf(struct_obj->fqresult_r2,"%s",PolyChar);
            } 
          }
          else{
            ksprintf(struct_obj->fqresult_r1,"%s",readadapt_rc_sam1);
            if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",readadapt_rc_sam2);} 
          }     
        }
      }
      else{
        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,seq_r2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          //fprintf(stderr,"TEST FQ BC%s\n",qual_r1);
          for(long unsigned int p = 0;p<strlen(seq_r1);p++){
            double dtemp1;double dtemp2;
            dtemp1 = mrand_pop(drand_alloc_nt);
            dtemp2 = mrand_pop(drand_alloc_nt);

            //fprintf(stderr,"NOT ADAPTER tmp1 %f\t%f\n",dtemp1,dtemp2);
            int base = seq_r1[p];
            int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
            qual_r1[p] = struct_obj->NtQual_r1[qscore];

            if (struct_obj->ErrorFlag == 'T'){
              fprintf(stderr,"ERROFRFLAG %c\n",struct_obj->ErrorFlag);
              double dtemp3;double dtemp4;
              dtemp3 = mrand_pop(drand_alloc_nt);
              dtemp4 = mrand_pop(drand_alloc_nt);
              //fprintf(stderr,"temp lol 3 %f \t temp 4 %f\n----------------\n",dtemp3,dtemp4);
              if (dtemp3 < struct_obj->NtErr_r1[qscore]){ErrorSub(dtemp4,seq_r1,p);}
            }
          }
          //fprintf(stderr,"TEST FQ AC%s\n",qual_r1);
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            for(int p = 0;p<seqlen;p++){
              double dtemp1;double dtemp2;
              dtemp1 = mrand_pop(drand_alloc_nt);
              dtemp2 = mrand_pop(drand_alloc_nt);

              int base = seq_r2[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);

              qual_r2[p] = struct_obj->NtQual_r2[qscore];
              if (struct_obj->ErrorFlag == 'T'){
                double dtemp3;double dtemp4;
                dtemp3 = mrand_pop(drand_alloc_nt);
                dtemp4 = mrand_pop(drand_alloc_nt);

                if (dtemp3 < struct_obj->NtErr_r2[qscore]){ErrorSub(dtemp4,seq_r2,p);}
              }
            }
            ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);}
        }
        if (struct_obj->SAMout){

          if(strcasecmp("true",struct_obj->QualFlag)==0){
            for(long unsigned int p = 0;p<strlen(seq_r1);p++){
              double dtemp1;double dtemp2;
              dtemp1 = mrand_pop(drand_alloc_nt);
              dtemp2 = mrand_pop(drand_alloc_nt);

              int base = seq_r1[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
              qual_r1[p] = struct_obj->NtQual_r1[qscore];

              if (struct_obj->ErrorFlag == 'T'){
                double dtemp3;double dtemp4;
                dtemp3 = mrand_pop(drand_alloc_nt);
                dtemp4 = mrand_pop(drand_alloc_nt);
                if (dtemp3 < struct_obj->NtErr_r1[qscore]){ErrorSub(dtemp4,seq_r1,p);}
              }
            }

            if (strcasecmp("PE",struct_obj->SeqType)==0){
              for(int p = 0;p<seqlen;p++){
                double dtemp1;double dtemp2;
                dtemp1 = mrand_pop(drand_alloc_nt);
                dtemp2 = mrand_pop(drand_alloc_nt);

                int base = seq_r2[p];

                int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
                qual_r2[p] = struct_obj->NtQual_r2[qscore];
                
                if (struct_obj->ErrorFlag == 'T'){
                  double dtemp3;double dtemp4;
                  dtemp3 = mrand_pop(drand_alloc_nt);
                  dtemp4 = mrand_pop(drand_alloc_nt);

                  if (dtemp3 < struct_obj->NtErr_r2[qscore]){ErrorSub(dtemp4,seq_r2,p);}
                }
              }
            }
          }
          //fprintf(stderr,"AFTER SEQ ERR \n%s \n",seq_r1);
          if (flag == 16 || flag == 81){DNA_complement(seq_r1);reverseChar(seq_r1,strlen(seq_r1));}
          else if (flag == 97){DNA_complement(seq_r2);reverseChar(seq_r2,strlen(seq_r2));}  
          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
          //fprintf(stderr,"SAVE OUTPUT \n%s \n",seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",seq_r2);}
        }
      }
      if (struct_obj->bgzf_fp1){
        if (struct_obj->fqresult_r1->l > 30000000){
          pthread_mutex_lock(&Fq_write_mutex);
          assert(bgzf_write(struct_obj->bgzf_fp1,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
          if (strcasecmp("PE",struct_obj->SeqType)==0){assert(bgzf_write(struct_obj->bgzf_fp2,struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->fqresult_r1->l =0;
          struct_obj->fqresult_r2->l =0;
        }
      }
      if (struct_obj->SAMout){
        if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
          n_cigar = 2;
          uint32_t cigar_bit_soft = bam_cigar_gen(strlen(struct_obj->fqresult_r1->s)-seqlen, BAM_CSOFT_CLIP);
          uint32_t cigar_arr[] = {cigar_bitstring,cigar_bit_soft};
          cigar = cigar_arr;

          if (strcasecmp("PE",struct_obj->SeqType)==0){
            uint32_t cigar_bit_soft2 = bam_cigar_gen(strlen(struct_obj->fqresult_r2->s)-seqlen, BAM_CSOFT_CLIP);
            uint32_t cigar_arr2[] = {cigar_bitstring,cigar_bit_soft2};
            cigar2 = cigar_arr2;
          }
        }
        else{
          n_cigar = 1;
          uint32_t cigar_arr[] = {cigar_bitstring};
          cigar = cigar_arr;
          cigar2 = cigar_arr;
        }
        
        size_t l_aux = 2; uint8_t mapq = 60;
        hts_pos_t min_beg, max_end, insert; //max_end, insert;
        min_beg = rand_start-struct_obj->size_cumm[chr_idx] - 1;
        max_end = rand_start-struct_obj->size_cumm[chr_idx] + fraglength - 1;
        insert = max_end - min_beg + 1;
        
        // since we change the orientation of the reads from the reverse strand, we also have to do it for the qual string
        // and it should only hold true for the sequence without the adapter
        if (flag == 16 || flag == 81){reverseChar(qual_r1,strlen(seq_r1));}
        else if (flag == 97){reverseChar(qual_r2,strlen(seq_r2));}  

        if (strcasecmp("PE",struct_obj->SeqType)==0){
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,
          n_cigar,cigar,chr_idx,max_end,insert,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag2,chr_idx,max_end,mapq,
          n_cigar,cigar2,chr_idx,min_beg,0-insert,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,qual_r2,l_aux);
        }
        else if (strcasecmp("SE",struct_obj->SeqType)==0){
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,
          n_cigar,cigar,-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
        }
        
        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) != 1);
          }
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->l = 0;
        }
        struct_obj->fqresult_r1->l =0;
        struct_obj->fqresult_r2->l =0;
      }
      
      //fprintf(stderr,"Thread %d \tlocal read number %d \t current read %d \t reads %d\n",struct_obj->threadno,localread,current_reads_atom,reads);
      //fprintf(stderr,"The seq %s \t qual %s \n",seq_r1,qual_r1);
      memset(qual_r1, 0, sizeof qual_r1); 
      memset(qual_r2, 0, sizeof qual_r2);  
      memset(seq_r1, 0, sizeof seq_r1);
      memset(seq_r1_mod, 0, sizeof seq_r1_mod);

      memset(seq_r2, 0, sizeof seq_r2);
      memset(seq_r2_mod, 0, sizeof seq_r2_mod);
      
      memset(readadapt, 0, sizeof readadapt);
      memset(readadapt2, 0, sizeof readadapt2);
            
      chr_idx = 0;
      iter++;
      localread++;
      current_reads_atom++;
    }
  }
  if (struct_obj->bgzf_fp1){
    if (struct_obj->fqresult_r1->l > 0){
      pthread_mutex_lock(&Fq_write_mutex);
      assert(bgzf_write(struct_obj->bgzf_fp1,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
      if (strcasecmp("PE",struct_obj->SeqType)==0){assert(bgzf_write(struct_obj->bgzf_fp2,struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
      pthread_mutex_unlock(&Fq_write_mutex);
      struct_obj->fqresult_r1->l =0;
      struct_obj->fqresult_r2->l =0;
    } 
  }

  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  for(int j=0; j<struct_obj->m;j++){bam_destroy1(struct_obj->list_of_reads[j]);}
  
  free(drand_alloc);
  free(drand_alloc_nt);
  free(drand_alloc_nt_adapt);
  free(drand_alloc_briggs);
  
  fprintf(stderr,"\t-> Number of reads generated by thread %d is %d \n",struct_obj->threadno,localread);

  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_flag,
                        const char* Sizefile,int FixedSize,int SizeDistType, int val1, int val2,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2, int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,const char* ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* MisMatchFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat,const char* Variant_flag,const char *VarType){
  //creating an array with the arguments to create multiple threads;
  //fprintf(stderr,"Random MacIntType %d\n",MacroRandType);

  int nthreads=thread_no;
  pthread_t mythreads[nthreads];

  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4; 

  int chr_total;
  char *genome_data;
  if (Specific_Chr[0] != NULL){
    while (Specific_Chr[chr_total]) {chr_total++;}
  }
  else{
    chr_total = faidx_nseq(seq_ref);
  }

  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_idx_arr[chr_total];
  size_t chr_size_cumm[chr_total+1];
  /*fprintf(stderr,"Chromosome count %d\n",chr_total);
  fprintf(stderr,"DONE WITH LOOP\n");*/
  
  if (chr_total < faidx_nseq(seq_ref)){
    for (int j = 0; j < faidx_nseq(seq_ref); j++){
      for (int i = 0; i < chr_total; i++){
        if(strcasecmp(faidx_iseq(seq_ref, j),Specific_Chr[i])==0){
          chr_idx_arr[i] = j;
          //fprintf(stderr,"index j %d and i %d and value %d\n",j,i,chr_idx_arr[i]);
        }
      } 
    }
  }
  else
  {
    for (int j = 0; j < faidx_nseq(seq_ref); j++){chr_idx_arr[j] = j;}
  }
  
  if(VCFformat != NULL && strcasecmp(Variant_flag,"bcf")==0){
    genome_data = full_vcf_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm,VCFformat,VarType);
  }
  else{
    if (chr_total == faidx_nseq(seq_ref)){
      genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
    }  
    else{
      genome_data = partial_genome_create(seq_ref,chr_total,chr_sizes,Specific_Chr,chr_size_cumm);
      for (int i = 0; i < chr_total; i++){
        chr_names[i] = Specific_Chr[i];
      }
    }
  }

  size_t genome_size = strlen(genome_data);;
  if (genome_data != NULL){
    fprintf(stderr,"\t-> Creating the large concatenated contig, with size of %lu bp\n",genome_size);
  
    Parsarg_for_Sampling_thread struct_for_threads[nthreads];

    // declare files and headers
    BGZF *bgzf_fp1 = NULL;
    BGZF *bgzf_fp2 = NULL;

    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader;
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));

    char file1[80];
    char file2[80];
    const char* fileprefix = OutputName; //"chr22_out";
    strcpy(file1,fileprefix);
    strcpy(file2,fileprefix);

    const char* suffix1;
    const char* suffix2;

    const char *mode;
    int alnformatflag = 0;
    if(strcasecmp("fa",OutputFormat)==0){
      //fprintf(stderr,"\t-> FA SE FILE\n");
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa";}
      else{suffix1 = "_r1.fa";suffix2 = "_r2.fa";}
    }
    else if(strcasecmp("fa.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa.gz";}
      else{suffix1 = "_r1.fa.gz";suffix2 = "_r2.fa.gz";}
    }
    else if(strcasecmp("fq",OutputFormat)==0){
      //fprintf(stderr,"\t-> FQ SE FILE\n");
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq";}
      else{suffix1 = "_r1.fq";suffix2 = "_r2.fq";}
    }
    else if(strcasecmp("fq.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq.gz";}
      else{suffix1 = "_r1.fq.gz";suffix2 = "_r2.fq.gz";}
    }
    else if(strcasecmp("sam",OutputFormat)==0){
      //fprintf(stderr,"\t-> SAM SE FILE\n");
      mode = "ws";
      suffix1 = ".sam";
      alnformatflag++;
    }
    else if(strcasecmp("bam",OutputFormat)==0){
      //fprintf(stderr,"\t-> BAM SE FILE\n");
      mode = "wb";//"wc";
      suffix1 = ".bam"; //".cram";
      alnformatflag++;
    }
    else if(strcasecmp("cram",OutputFormat)==0){
      mode = "wc";
      suffix1 = ".cram";
      alnformatflag++;
    }
    else{fprintf(stderr,"\t-> Fileformat is currently not supported \n");}
    strcat(file1,suffix1);

    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(alnformatflag == 0){
      //fprintf(stderr,"not bam loop \n");
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp1 = bgzf_open(filename1,mode);
      bgzf_mt(bgzf_fp1,mt_cores,bgzf_buf);
      
      //fprintf(stderr,"\t-> BGZF FILE\n");
      if(strcasecmp("PE",SeqType)==0){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp2 = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp2,mt_cores,bgzf_buf);
      }
    }
    else{
      fprintf(stderr,"Fasta input file name is %s\n",FastaFileName);
      char *ref =(char*) malloc(10 + strlen(FastaFileName) + 1);
      sprintf(ref, "reference=%s", FastaFileName);
      //char *ref =(char*) malloc(10 + strlen("Test_Examples/Mycobacterium_leprae.fa.gz") + 1);
      //sprintf(ref, "reference=%s", "Test_Examples/Mycobacterium_leprae.fa.gz");
      hts_opt_add((hts_opt **)&fmt_hts->specific,ref);
      fprintf(stderr,"Writing mode is %s\n",mode);
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();
      
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,seq_ref,chr_total,chr_idx_arr,genome_size);
    }
    //fprintf(stderr,"\t-> AFTER OUTPUT FORMAT\n");

    int number; int* Frag_len; double* Frag_freq;
    fprintf(stderr,"SizeDistType %d \t FixedSize %d\n",SizeDistType,FixedSize);
    if(FixedSize==-1 && SizeDistType==-1){
      fprintf(stderr,"\t-> FRAG DIST FILE\n");
      Frag_len = new int[LENS];
      Frag_freq = new double[LENS];
      FragArray(number,Frag_len,Frag_freq,Sizefile); //Size_dist_sampling //"Size_dist/Size_freq_modern.txt"
      //fprintf(stderr,"\t-> FRAG ARRAY LE\n");
    }
    else if(SizeDistType!=-1){
      fprintf(stderr,"\t-> FRAG DIST LENGTH\n");
      Frag_len = new int[LENS];
      Frag_freq = new double[LENS];
      FragDistArray(number,Frag_len,Frag_freq,SizeDistType,seed,val1, val2);
    }
    else if(FixedSize!=-1){number = -1;}
    fprintf(stderr,"THE NUMBER IS %d\n",number);
    const char *freqfile_r1; //"Qual_profiles/AccFreqL150R1.txt";
    const char *freqfile_r2;
    int outputoffset = qualstringoffset;
    unsigned long readcyclelength;
    //fprintf(stderr,"\t-> FRAG ARRAY LE\n");
    ransampl_ws ***QualDist;
    char nt_qual_r1[1024];
    ransampl_ws ***QualDist2;
    char nt_qual_r2[1024];
    //fprintf(stderr,"\t-> QUAL POINTER POINTER POINTER\n");
    double ErrArray_r1[1024];
    double ErrArray_r2[1024];
    //fprintf(stderr,"\t-> BEFORE QUAL STRING IF\n");
    if(strcasecmp("true",QualStringFlag)==0){ //|| strcasecmp("bam",OutputFormat)==0
      freqfile_r1 = QualProfile1;
      QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1,readcyclelength);
      //fprintf(stderr,"\t-> CREATING QUALDIST\n");
      if(strcasecmp("PE",SeqType)==0){
        //fprintf(stderr,"\t-> PE LOOP\n");
        freqfile_r2 = QualProfile2;
        QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2,readcyclelength);
      }
    }
    //fprintf(stderr,"\t-> AFTER QUAL STRING IF\n");
    
    int maxsize = 20;
    char polynucleotide;
    //fprintf(stderr,"\t-> BEFORE POLY\n");
    if (Polynt != NULL && strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}
    //fprintf(stderr,"\t-> AFTER POLY\n");

    const char *Sub_mat = SubProfile;
    double* DeamFreqArray = new double[LENS];
    int deamcyclelength = 0;
    if (SubProfile != NULL){
      std::cout << Sub_mat << std::endl;
      DeamFreqArray = DeamFileArray(DeamFreqArray,SubProfile,deamcyclelength);
      std::cout << DeamFreqArray[0] << std::endl;
      fprintf(stderr,"INFERRED READ LENGTH %d \n",deamcyclelength);
    }
    //else{std::cout << "LOLRT " << std::endl;}
    
    
    //fprintf(stderr,"\t-> Before threads\n");
    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;

      struct_for_threads[i].fqresult_r2 =new kstring_t;
      struct_for_threads[i].fqresult_r2 -> l = 0;
      struct_for_threads[i].fqresult_r2 -> m = 0;
      struct_for_threads[i].fqresult_r2 -> s = NULL;

      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_idx_array = chr_idx_arr;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].RandMacro = RandMacro;

      struct_for_threads[i].FragLen = Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = number;
      struct_for_threads[i].FixedSize = FixedSize;

      struct_for_threads[i].NtQual_r1 = nt_qual_r1;
      struct_for_threads[i].NtQual_r2 = nt_qual_r2;
      struct_for_threads[i].QualDist_r1 = QualDist;
      struct_for_threads[i].QualDist_r2 = QualDist2;
      struct_for_threads[i].NtErr_r1 = ErrArray_r1;
      struct_for_threads[i].NtErr_r2 = ErrArray_r2;

      struct_for_threads[i].MisMatch = DeamFreqArray;
      struct_for_threads[i].SubFlag = MisMatchFlag;
      struct_for_threads[i].MisLength = (int) deamcyclelength;
      struct_for_threads[i].readcycle = (int) readcyclelength;
      struct_for_threads[i].reads = reads;

      struct_for_threads[i].bgzf_fp1 = bgzf_fp1;
      struct_for_threads[i].bgzf_fp2 = bgzf_fp2;
      struct_for_threads[i].SAMout = SAMout;
      struct_for_threads[i].SAMHeader = SAMHeader;
      struct_for_threads[i].l = 0;
      struct_for_threads[i].m = maxsize;
      struct_for_threads[i].list_of_reads = (bam1_t**) malloc(sizeof(bam1_t)*maxsize); // need to free this space

      for(int j=0; j<maxsize;j++){struct_for_threads[i].list_of_reads[j]=bam_init1();} // but also destroy the bam_init1 objects    

      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].Adapter_2 = Adapter_2;
      struct_for_threads[i].Briggs_flag = Briggs_flag;
      struct_for_threads[i].BriggsParam = BriggsParam;
      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].SeqType = SeqType;
      struct_for_threads[i].QualFlag = QualStringFlag;
      struct_for_threads[i].PolyNt = polynucleotide;
      struct_for_threads[i].ErrorFlag = (char) ErrorFlag[0];

      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (size_t*)malloc(sizeof(size_t) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Sampling_threads,&struct_for_threads[i]);
    }
    
    for (int i = 0; i < nthreads; i++){  
      pthread_join(mythreads[i],NULL);
    }
        
    if(alnformatflag == 0){
      bgzf_close(bgzf_fp1);
      if(strcasecmp("PE",SeqType)==0){bgzf_close(bgzf_fp2);}
    }
    else{
      sam_hdr_destroy(SAMHeader);
      sam_close(SAMout);
    } 
    
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);
      free(struct_for_threads[i].list_of_reads);
      delete struct_for_threads[i].fqresult_r1;

      free(struct_for_threads[i].fqresult_r2 -> s);
      delete struct_for_threads[i].fqresult_r2;      
    }
    
    if(strcasecmp("true",QualStringFlag)==0){
      for(int b=0;b<5;b++){
        for(int pos = 0 ; pos< (int) readcyclelength;pos++){
          ransampl_free(QualDist[b][pos]);
        }
        delete[] QualDist[b];
      }
      delete[] QualDist;

      if(strcasecmp("PE",SeqType)==0){
        for(int b=0;b<5;b++){
          for(int pos = 0 ; pos< (int) readcyclelength;pos++){
            ransampl_free(QualDist2[b][pos]);
          }
          delete[] QualDist2[b];
        }
        delete[] QualDist2;
      }
    }
    
    free(fmt_hts);
    if(Sizefile!=NULL){
      delete[] Frag_freq;
      delete[] Frag_len;
    }
    
    if(SubProfile != NULL){
      delete[] DeamFreqArray;
    }
    
    free(genome_data);
    fflush(stderr);
  }
  return NULL;
}
