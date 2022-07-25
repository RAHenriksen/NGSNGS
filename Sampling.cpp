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
#include <htslib/thread_pool.h>

#include <pthread.h>

#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "NGSNGS_func.h"
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "ThreadGeneration.h"
#include "Sampling.h"

#define LENS 4096
#define MAXBINS 100
unsigned char nuc2int[255];

pthread_mutex_t Fq_write_mutex;

int MacroRandType = 0;

void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  int MacroRandType;
  MacroRandType = struct_obj -> RandMacro;
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  size_t genome_len = strlen(struct_obj->genome);
  mrand_t *drand_alloc = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_nt = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_nt_adapt = mrand_alloc(MacroRandType,loc_seed);
  mrand_t *drand_alloc_briggs = mrand_alloc(MacroRandType,loc_seed);
   
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

  size_t reads = struct_obj -> reads;
  size_t BufferLength = struct_obj -> BufferLength;

  size_t rand_start;

  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};

  size_t localread = 0;
  int iter = 0;
  size_t current_reads_atom = 0;
  int readsizelimit;

  extern int SIG_COND;
  while (current_reads_atom < reads &&SIG_COND) {
    double rand_val = mrand_pop(drand_alloc);
    rand_start = rand_val * (genome_len-300)+1; //genome_len-100000;
    double rand_val_len = mrand_pop(drand_alloc);
    //std::cout << rand_val_len << std::endl;
    // Fragment length creation
    int fraglength;
    if (struct_obj->No_Len_Val != -1){
      // random start and length are not dependent on the same rand val
      //fprintf(stderr,"Legth%f \n",mrand_pop(drand_alloc));
      int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val_len);
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
    if (fraglength > readsizelimit){strncpy(seq_r1,struct_obj->genome+rand_start-1,readsizelimit);}   // case 1
    else {strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    
    if(strcasecmp("PE",struct_obj->SeqType)==0){
      if (fraglength > readsizelimit){

        strncpy(seq_r2,struct_obj->genome+rand_start+fraglength-1-readsizelimit,readsizelimit);
        } // case 1
      else {strncpy(seq_r2,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    }

    //fprintf(stderr,"Chrindx %d \t chromosome name %s \t chromosome length %zu \t cumulative length %zu \t start pos %zu \n",chr_idx,struct_obj->names[chr_idx],struct_obj->size_cumm[chr_idx+1]-struct_obj->size_cumm[chr_idx],struct_obj->size_cumm[chr_idx+1],rand_start);

    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N'); //first encounter with N
    pch2 = strrchr(seq_r1,'N'); //last encounter with N

    int seqlen = strlen(seq_r1);
    size_t n_cigar;const uint32_t *cigar;const uint32_t *cigar2; uint16_t flag; uint16_t flag2;
    uint32_t cigar_bitstring = bam_cigar_gen(seqlen, BAM_CMATCH);

    // for unaligned 
    size_t n_cigar_unmap=1;const uint32_t *cigar_unmap;
    uint32_t cigar_bitstring_unmap = bam_cigar_gen(seqlen, BAM_CSOFT_CLIP);
    uint32_t cigar_arr_unmap[] = {cigar_bitstring_unmap};
    cigar_unmap = cigar_arr_unmap;

    if ((int )(pch-seq_r1+1) == 1 || (int)(pch2-seq_r1+1)  == seqlen){
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);}
    else if ((rand_start-struct_obj->size_cumm[chr_idx+1])<fraglength){
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);
    }
    else if (rand_start+fraglength>struct_obj->size_cumm[chr_idx+1]){
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);
    }
    else{
      // then all the same start pos would have the same strand no matter the chromosome?
      int strand = (int) (rand_start%2);//(int) rand_r(&loc_seed)%2;//1;//rand() % 2;
      
      //fprintf(stderr,"STRAND %d \t random start %zu \t random start int %d \t modulo %zu\n",strand,rand_start,(int) rand_start, rand_start%2);

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
          MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
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
          MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          // IS THIS THE WAY TO DO IT?
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          }
          //fprintf(stderr,"SEQUENCE %s\n",seq_r1);
        }
      }
      
      char READ_ID[1024]; int read_id_length;
    
      char *dummydummy = struct_obj->names[chr_idx];
      if (struct_obj->Variant_flag!=NULL && strcasecmp(struct_obj->Variant_flag ,"bcf")==0)
	    dummydummy = struct_obj->names[0];
      
      read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%zu-%zu_length:%d", struct_obj->threadno, rand_id,strand,dummydummy,
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
          ksprintf(struct_obj->fqresult_r1,">%s_R1\n%s\n",READ_ID,readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s_R2\n%s\n",READ_ID,readadapt2);}
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
            ksprintf(struct_obj->fqresult_r1,"@%s_R1\n%s\n+\n%s\n",READ_ID,PolyChar,QualPoly);
          }
          else{ksprintf(struct_obj->fqresult_r1,"@%s_R1\n%s\n+\n%s\n",READ_ID,readadapt,qual_r1);}
          
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
              ksprintf(struct_obj->fqresult_r2,"@%s_R1\n%s\n+\n%s\n",READ_ID,PolyChar,QualPoly);
            }
            else{ksprintf(struct_obj->fqresult_r2,"@%s_R2\n%s\n+\n%s\n",READ_ID,readadapt2,qual_r2);}
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
            sprintf(readadapt_err1, "%*s", (int)strlen(readadapt) - (int)strlen(seq_r1), readadapt+(int)strlen(seq_r1));
            sprintf(read_rc_sam1, "%.*s", (int)strlen(seq_r1), readadapt);
            DNA_complement(read_rc_sam1);reverseChar(read_rc_sam1,strlen(read_rc_sam1));
            //fprintf(stderr,"string before \t\t %s \n",readadapt_rc_sam1);
            sprintf(readadapt_rc_sam1, "%s", read_rc_sam1);
            //fprintf(stderr,"string after \t\t %s \n",readadapt_rc_sam1);
            strcat(readadapt_rc_sam1, readadapt_err1);
            //fprintf(stderr,"string after II \t %s \n ----------------- \n",readadapt_rc_sam1);
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              // the mate to flag 81 will be the second in pair on the forward strand - just to esnrue identical names
              sprintf(readadapt_rc_sam2, "%s",readadapt2);
            } 
          }
          else if (flag == 97){
            // 97 first in pair forward strand only for Paired end.
            sprintf(readadapt_rc_sam1, "%s",readadapt); //to ensure all the reads have identifcal names for the sam output

            // the mate is reverse
            sprintf(readadapt_err2, "%*s", (int)strlen(readadapt2)-(int)strlen(seq_r2), readadapt2+(int)strlen(seq_r2));
            sprintf(read_rc_sam2, "%.*s", (int)strlen(seq_r2), readadapt2); // copy the sequence 2 from the sequence 2 + adapter
            // type cast to int
            DNA_complement(read_rc_sam2);reverseChar(read_rc_sam2,strlen(read_rc_sam1)); // reverse complement sequence 2 
            
            //fprintf(stderr,"string before \t\t %s \n",readadapt_rc_sam2);
            sprintf(readadapt_rc_sam2, "%s", read_rc_sam2);
            //fprintf(stderr,"string after I \t\t %s \n",readadapt_rc_sam2);
            strcat(readadapt_rc_sam2, readadapt_err2);
            //fprintf(stderr,"string after II\t\t %s \n",readadapt_rc_sam2);
            //sprintf(readadapt_rc_sam2, "%s%s", read_rc_sam2,readadapt_err2); // create rev comp sequence 2 + adapter

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
          ksprintf(struct_obj->fqresult_r1,">%s_R1\n%s\n",READ_ID,seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s_R2\n%s\n",READ_ID,seq_r2);}
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
              //fprintf(stderr,"ERROFRFLAG %c\n",struct_obj->ErrorFlag);
              double dtemp3;double dtemp4;
              dtemp3 = mrand_pop(drand_alloc_nt);
              dtemp4 = mrand_pop(drand_alloc_nt);
              //fprintf(stderr,"temp lol 3 %f \t temp 4 %f\n----------------\n",dtemp3,dtemp4);
              if (dtemp3 < struct_obj->NtErr_r1[qscore]){ErrorSub(dtemp4,seq_r1,p);}
            }
          }
          //fprintf(stderr,"TEST FQ AC%s\n",qual_r1);
          ksprintf(struct_obj->fqresult_r1,"@%s_R1\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
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
            ksprintf(struct_obj->fqresult_r2,"@%s_R2\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);}
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
        if (struct_obj->fqresult_r1->l > BufferLength){
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

        char READIDR1[1024];char READIDR2[1024];
        const char* suffR1 = "_R1";const char* suffR2 = "_R2";
        strcpy(READIDR1,READ_ID);strcat(READIDR1,suffR1);
        if (strcasecmp("PE",struct_obj->SeqType)==0){
          strcpy(READIDR2,READ_ID);strcat(READIDR2,suffR2);
          if (struct_obj->NoAlign == 'T'){
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR1),READIDR1,4,-1,-1,
            255,n_cigar_unmap,cigar_unmap,-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,NULL,0);
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR2),READIDR2,4,-1,-1,
            255,n_cigar_unmap,cigar_unmap,-1,-1,0,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,NULL,0);
          }
          else{
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR1),READIDR1,flag,chr_idx,min_beg,mapq,
            n_cigar,cigar,chr_idx,max_end,insert,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR2),READIDR2,flag2,chr_idx,max_end,mapq,
            n_cigar,cigar2,chr_idx,min_beg,0-insert,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,qual_r2,l_aux);
          }
        }
        else if (strcasecmp("SE",struct_obj->SeqType)==0){
          //fprintf(stderr,"READID v2 %s\n",READIDR1);
          if (struct_obj->NoAlign == 'T'){
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR1),READIDR1,4,-1,-1,
            255,n_cigar_unmap,cigar_unmap,-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,NULL,0);
          }
          else{
            bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length+strlen(suffR1),READIDR1,flag,chr_idx,min_beg,mapq,
            n_cigar,cigar,-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
          }
          //const char* MDtag = "\tMD:Z:";
          //bam_aux_update_str(struct_obj->list_of_reads[struct_obj->l-1],"MD",5,"MDtag");
        }

        int assert_int = 0;
        if(strcasecmp("cram",struct_obj->OutputFormat)==0){assert_int = 1;}

        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            //fprintf(stderr,"INSIDE FOR LOOP WITH ASSERTION\n");
            assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) != assert_int);
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
  
  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);

  pthread_exit(NULL);
}
