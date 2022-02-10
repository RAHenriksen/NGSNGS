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

#define LENS 4096
#define MAXBINS 100
char nuc2int[255];

pthread_mutex_t Fq_write_mutex;

// ---------------------- SINGLE-END ---------------------- //
struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  size_t *size_cumm;
  const char **names;

  int* FragLen;
  double* FragFreq;
  int No_Len_Val;

  char *NtQual_r1;
  char *NtQual_r2;
  ransampl_ws ***QualDist_r1; //double* Qualfreq;
  ransampl_ws ***QualDist_r2; //double* Qualfreq;
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

  char PolyNt;
};
      
void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;
  //fprintf(stderr,"%s",struct_obj->SeqType);
  // creating random objects for all distributions.
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  size_t genome_len = strlen(struct_obj->genome);
  
  //srand48(loc_seed);
  struct drand48_data buffer;
  srand48_r(loc_seed, &buffer);
  /*for (size_t i = 0; i < 100; i++){
    double dtemp1;double dtemp2;
    drand48_r(&buffer, &dtemp1);
    drand48_r(&buffer, &dtemp2);
    fprintf(stderr,"Double value %f \t %f \n ----------\n", dtemp1,dtemp2);
  }
  exit(0);*/
  
  //coverage method2
  char seq_r1[1024] = {0}; //original
  char seq_r1_mod[1024] = {0}; //modifying


  char read[1024] = {0};
  char readadapt[1024] = {0};

  char seq_r2[1024] = {0};
  char seq_r2_mod[1024] = {0};
  char read2[1024] = {0};
  char readadapt2[1024] = {0};
  
  int reads = struct_obj -> reads;

  //float cov_current = 0;
  size_t rand_start;
  //int nread = 0;

  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};

  int localread = 0;
  int iter = 0;
  int current_reads_atom = 0;
  int readsizelimit;

  double ErrArray[8] = {0.6309573,0.2511886,0.03162278,0.006309573,0.001995262,0.0005011872,0.0001995262,0.0001};

  double dtemp1;double dtemp2;
  
  int C_count = 0;int G_count = 0;int A_count = 0;int T_count = 0;double CT_count = 0;double GA_count = 0;
  
  while (current_reads_atom < reads){
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX); //random between 0 and 1
    double rand_val2 = rand_val*RAND_MAX; // between 0 and maximum 
    double rand_val3 = myrand((unsigned int) rand_val2); //((double) rand_r(&test)/ RAND_MAX);// 
    rand_start = rand_val3 * (genome_len-300); //genome_len-100000;
    // fprintf(stderr,"READ %d with rand val 1: %f 2: %lf 3: %f\n",current_reads_atom,rand_val,rand_val2,rand_val3); // ALL RANDOM
    // Fragment length creation
    int fraglength;
    if (struct_obj->No_Len_Val != -1){
      int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val);
      fraglength =  struct_obj->FragLen[lengthbin]; //struct_obj->FragLen[lengthbin];//75;
    }
    else{
      //fprintf(stderr,"FIXED LENGTH \n");
      fraglength = struct_obj->FixedSize;
    } 
    if(strcasecmp("false",struct_obj->QualFlag)==0){
      readsizelimit = fraglength;
      // fprintf(stderr,"READ SIZE LIMIT IF %d \n",readsizelimit);
    }
    else{
      readsizelimit = struct_obj->readcycle;
      // fprintf(stderr,"READ SIZE LIMIT ELSE %d\n",readsizelimit);
    }
    
    //fprintf(stderr,"Fragment length %d\n",fraglength);

    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;} //extracting the correct chromosome id for the start pos

    if (fraglength > readsizelimit){strncpy(seq_r1,struct_obj->genome+rand_start-1,readsizelimit);}   // case 1
    else {strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    
    if(strcasecmp("PE",struct_obj->SeqType)==0){
      if (fraglength > readsizelimit){
        /*std::cout << rand_start+fraglength-1-readsizelimit << std::endl;
        std::cout << rand_start-1 << std::endl;*/
        strncpy(seq_r2,struct_obj->genome+rand_start+fraglength-1-readsizelimit,readsizelimit);
        } // case 1
      else {strncpy(seq_r2,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    }

    int rand_id = (rand_val * fraglength-1); //100
    double rand_val4 = myrand(rand_id+1);
    double rand_val5 = myrand(rand_id+2);
    //fprintf(stderr,"READ ID %d \t with rand val 4: %f 5: %f\n",rand_id,rand_val4,rand_val5); //All random

    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N');
    pch2 = strrchr(seq_r1,'N');
    //fprintf(stderr,"STRCHR %d %s\n",(int)(pch2-seq_r1+1),seq_r1);
    //if (pch != NULL){continue;}
    int seqlen = strlen(seq_r1);
    int seqlen2 = strlen(seq_r2);
    
    size_t n_cigar;const uint32_t *cigar;const uint32_t *cigar2; uint16_t flag; uint16_t flag2;
    uint32_t cigar_bitstring = bam_cigar_gen(seqlen, BAM_CMATCH);
    //fprintf(stderr,"SEQUENCE LENGTH %d \t %d \t %d\n", strlen(seq_r1),seqlen,strlen(seq_r2));
    if ((int )(pch-seq_r1+1) == 1 || (int)(pch2-seq_r1+1)  == seqlen){
      memset(seq_r1, 0, sizeof seq_r1);memset(seq_r2, 0, sizeof seq_r2);}
    else{
      int strand = (int) rand_r(&loc_seed)%2;//1;//rand() % 2;
      //fprintf(stderr,"strand %d \n",strand); // 1 or 0

      if (struct_obj->SAMout){
        // in bam files the reads are all aligning to the forward strand, but the flags identify read property
        if (strcasecmp("SE",struct_obj->SeqType)==0){
          if (strand == 0){flag = 0;} // se read forward strand
          else{flag = 16;
            DNA_complement(seq_r1);
            reverseChar(seq_r1);
          } // reverse strand
        }
        if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 0){
          flag = 97;  //Paired, mate reverse, first
          flag2 = 145; // Paired, reverse strand, second
          DNA_complement(seq_r2);
          reverseChar(seq_r2);
        }
        else if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 1){
          flag = 81; // Paired, reverse strand, first
          flag2 = 161; // Paired, mate reverse, second
          DNA_complement(seq_r1);
          reverseChar(seq_r1);
        }

        if(strcasecmp(struct_obj->Briggs_flag,"true")==0){
          // SimBriggsModel(seq_r1, seq_r1_mod, fraglength, 0.024, 0.36, 0.68, 0.0097,loc_seed);
          //if (seq_r1[0]=='C' ||seq_r1[0]=='c'){C_count++;} // all C
          //if (seq_r1[strlen(seq_r1)-1]=='G'||seq_r1[strlen(seq_r1)-1]=='g'){G_count++;}

          SimBriggsModel(seq_r1, seq_r1_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,buffer);

          //if ((seq_r1[0]=='C'||seq_r1[0]=='c') && (seq_r1_mod[0]=='T' || seq_r1_mod[0]=='t')){CT_count++;}
          //if ((seq_r1[strlen(seq_r1)-1]=='G'||seq_r1[strlen(seq_r1)-1]=='g') && (seq_r1_mod[strlen(seq_r1)-1]=='A' || seq_r1_mod[strlen(seq_r1)-1]=='a')){GA_count++;}
          
          strncpy(seq_r1, seq_r1_mod, sizeof(seq_r1));
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            SimBriggsModel(seq_r2, seq_r2_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,buffer);
            strncpy(seq_r2, seq_r2_mod, sizeof(seq_r2));
          }
        }

        if (flag == 16 || flag == 81){DNA_complement(seq_r1);reverseChar(seq_r1);}
        else if (flag == 97){DNA_complement(seq_r2);reverseChar(seq_r2);}  
      }
      else{
        // in fasta and fastq the sequences need to be on forward or reverse strand, i.e we need reverse complementary        
        if (strcasecmp("SE",struct_obj->SeqType)==0 && strand == 1){
          DNA_complement(seq_r1);
          reverseChar(seq_r1);
          //if (seq_r1[strlen(seq_r1)-1]=='G'){G_count++;}
        }

        if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 0){
          DNA_complement(seq_r2);
          reverseChar(seq_r2);
        }
        else if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 1){
          DNA_complement(seq_r1);
          reverseChar(seq_r1);
        }

        // ADDING PMD AFTER THE STRAND IS CHOSEN AS TO NOT INFLUENCE THE SYMMETRY MORE
        if(strcasecmp(struct_obj->Briggs_flag,"true")==0){
          // SimBriggsModel(seq_r1, seq_r1_mod, fraglength, 0.024, 0.36, 0.68, 0.0097,loc_seed);
          SimBriggsModel(seq_r1, seq_r1_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,buffer);
          strncpy(seq_r1, seq_r1_mod, sizeof(seq_r1));
          
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            SimBriggsModel(seq_r2, seq_r2_mod, fraglength,struct_obj->BriggsParam[0], 
                                                        struct_obj->BriggsParam[1], 
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,buffer);
            strncpy(seq_r2, seq_r2_mod, sizeof(seq_r2));
          }
        }
      }
      
      char READ_ID[1024]; int read_id_length;
      read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%d-%d_length:%d", struct_obj->threadno, rand_id,strand,
        struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
        fraglength);
      
      /*read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%d_length:%d", struct_obj->threadno, rand_id*2,strand,
        rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],fraglength);*/
      
      if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
        //consider stpcpy for concatenation
        strcpy(read, seq_r1);
        strcat(read,struct_obj->Adapter_1);
        //fprintf(stderr,"STRING LENGTH SEQ %d AND ADAPTER %d AND READ %d\n",strlen(seq_r1),strlen(struct_obj->Adapter_1),strlen(read));
        //fprintf(stderr,"READ %s\n",read);
        if(strcasecmp("false",struct_obj->QualFlag)==0){
          //fprintf(stderr,"INSIDE false\n");
          strncpy(readadapt, read, strlen(read));
          readsizelimit = strlen(read);
        }
        else{strncpy(readadapt, read, readsizelimit);}
        //fprintf(stderr,"INSIDE ADAPTER TRUE\n");
        //fprintf(stderr,"READ SIZE LIMIT %d\n",readsizelimit);

        if (strcasecmp("PE",struct_obj->SeqType)==0){
          strcpy(read2, seq_r2);
          strcat(read2,struct_obj->Adapter_2);
          if(strcasecmp("false",struct_obj->QualFlag)==0){
            //fprintf(stderr,"INSIDE false\n");
            strncpy(readadapt2, read2, strlen(read2));
            readsizelimit = strlen(read2);
          }
          else{strncpy(readadapt2, read2, readsizelimit);}
          //fprintf(stderr,"INSIDE PE TRUE\n");
        }

        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          //fprintf(stderr,"INSIDE FA TRUE\n");
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,readadapt2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          //fprintf(stderr,"INSIDE FQ TRUE\n");
          for(int p = 0;p<strlen(readadapt);p++){
            double dtemp1;double dtemp2;
            drand48_r(&buffer, &dtemp1);
            drand48_r(&buffer, &dtemp2);
            
            int base = seq_r1[p];
            int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
            qual_r1[p] = struct_obj->NtQual_r1[qscore];
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
            for(int p = 0;p<strlen(readadapt2);p++){
              double dtemp1;double dtemp2;
              drand48_r(&buffer, &dtemp1);
              drand48_r(&buffer, &dtemp2);

              int base = seq_r2[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
              qual_r2[p] = struct_obj->NtQual_r2[qscore];
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
          //fprintf(stderr,"INSIDE SAM TRUE\n");
          if(strcasecmp("true",struct_obj->QualFlag)==0){
            //fprintf(stderr,"INSIDE QUAL STRING IF\n");
            for(int p = 0;p<strlen(readadapt);p++){
              double dtemp1;double dtemp2;
              drand48_r(&buffer, &dtemp1);
              drand48_r(&buffer, &dtemp2);
              int base = seq_r1[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);
              qual_r1[p] = struct_obj->NtQual_r1[qscore];
            }
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              //fprintf(stderr,"INSIDE SAM PE TRUE\n");
              for(int p = 0;p<strlen(readadapt2);p++){
                double dtemp1;double dtemp2;
                drand48_r(&buffer, &dtemp1);
                drand48_r(&buffer, &dtemp2);
                int base = seq_r2[p];
                int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
                qual_r2[p] = struct_obj->NtQual_r2[qscore];
              }
            }
          }
          ksprintf(struct_obj->fqresult_r1,"%s",readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",readadapt2);}
        }
      }
      else{
        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,seq_r2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          for(int p = 0;p<seqlen;p++){
            drand48_r(&buffer, &dtemp1);
            drand48_r(&buffer, &dtemp2);
            //fprintf(stderr,"Read %d SE rand val 4: %f and rand val 5: %f\n",current_reads_atom,myrand((unsigned int) (rand_val-p-current_reads_atom)),myrand((unsigned int) (rand_val+p+current_reads_atom)));
            int base = seq_r1[p];
            int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);

            //int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],myrand((unsigned int) (rand_val-p-current_reads_atom)),myrand((unsigned int) (rand_val+p+current_reads_atom)));
            qual_r1[p] = struct_obj->NtQual_r1[qscore];

            drand48_r(&buffer, &dtemp1);
            drand48_r(&buffer, &dtemp2);
            /*if (dtemp1 < struct_obj->NtErr_r1[qscore]){
              if (dtemp2 <= 0.25){seq_r1[p] = 'X';}
              else if (0.25 < dtemp2 && dtemp2 <= 0.5){seq_r1[p] = 'Q';}
              else if (0.5 < dtemp2 && dtemp2 <= 0.75){seq_r1[p] = 'Z';}
              else if (0.75 < dtemp2 && dtemp2 <= 1){seq_r1[p] = 'W';}
            }*/
            
            //std::cout << myrand((unsigned int) (rand_val-p-current_reads_atom)) << std::endl;
            //fprintf(stderr,"Double value SE %f \t %f \n", dtemp1,dtemp2);
            //std::cout << "-----------------------------" << std::endl;

          }
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            for(int p = 0;p<seqlen;p++){
              drand48_r(&buffer, &dtemp1);
              drand48_r(&buffer, &dtemp2);
              int base = seq_r2[p];
              int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
              //int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],myrand((unsigned int) (rand_val2+p+current_reads_atom)),myrand((unsigned int) (rand_val2-p-current_reads_atom)));
              qual_r2[p] = struct_obj->NtQual_r2[qscore];
              drand48_r(&buffer, &dtemp1);
              drand48_r(&buffer, &dtemp2);
              /*if (dtemp1 < struct_obj->NtErr_r2[qscore]){
                if (dtemp2 <= 0.25){seq_r2[p] = 'Y';}
                else if (0.25 < dtemp2 && dtemp2 <= 0.5){seq_r2[p] = 'R';}
                else if (0.5 < dtemp2 && dtemp2 <= 0.75){seq_r2[p] = 'O';}
                else if (0.75 < dtemp2 && dtemp2 <= 1){seq_r2[p] = 'P';}
              }*/
            }
            ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);}
        }
        if (struct_obj->SAMout){
          //fprintf(stderr,"SAM LOOP \n");
          if(strcasecmp("true",struct_obj->QualFlag)==0){
            //fprintf(stderr,"QUAL DIST IF \n");
            for(int p = 0;p<seqlen;p++){
              double dtemp1;double dtemp2;
              drand48_r(&buffer, &dtemp1);
              drand48_r(&buffer, &dtemp2);
              int base = seq_r1[p];

              int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],dtemp1,dtemp2);              
              //int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],myrand((unsigned int) (rand_val-p-current_reads_atom)),myrand((unsigned int) (rand_val+p+current_reads_atom)));
              qual_r1[p] = struct_obj->NtQual_r1[qscore];
            }
            // Read_Qual_new(seq_r1,qual,loc_seed,struct_obj->Qualfreq,0);
            if (strcasecmp("PE",struct_obj->SeqType)==0){
              for(int p = 0;p<seqlen;p++){
                double dtemp1;double dtemp2;
                drand48_r(&buffer, &dtemp1);
                drand48_r(&buffer, &dtemp2);
                int base = seq_r2[p];

                //int qscore = ransampl_draw2(struct_obj->QualDist_r1[nuc2int[base]][p],myrand((unsigned int) (rand_val-p-current_reads_atom)),myrand((unsigned int) (rand_val+p+current_reads_atom)));
                int qscore = ransampl_draw2(struct_obj->QualDist_r2[nuc2int[base]][p],dtemp1,dtemp2);
                qual_r2[p] = struct_obj->NtQual_r2[qscore];
              }
              // Read_Qual_new(seq_r2,qual2,loc_seed,struct_obj->Qualfreq,0);
              
            }
          }
          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",seq_r2);}
        }
      }
      if (struct_obj->bgzf_fp1){
        if (struct_obj->fqresult_r1->l > 30000000){
          //fprintf(stderr,"\t Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
          pthread_mutex_lock(&Fq_write_mutex);
          // bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
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
        if (strcasecmp("PE",struct_obj->SeqType)==0){
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,
          n_cigar,cigar,chr_idx,max_end,insert,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag2,chr_idx,max_end,mapq,
          n_cigar,cigar2,chr_idx,min_beg,0-insert,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,qual_r2,l_aux);
        }
        else if (strcasecmp("SE",struct_obj->SeqType)==0){
          bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,
          n_cigar,cigar,-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
          
          /*uint8_t *old_md = bam_aux_get(struct_obj->list_of_reads[struct_obj->l++], "MD");
          bam_aux_del(struct_obj->list_of_reads[struct_obj->l++], old_md);
          bam_aux_append(struct_obj->list_of_reads[struct_obj->l++], "MD", 'Z', struct_obj->fqresult_r1->l + 1, (uint8_t*)struct_obj->fqresult_r1->s);
          if (bam_aux_append(struct_obj->list_of_reads[struct_obj->l++], "MD", 'Z', struct_obj->fqresult_r1->l + 1, (uint8_t*)struct_obj->fqresult_r1->s) < 0)
            fprintf(stderr,"AUXILLARY FAILED\n");*/
        }
        
        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) != 1);
            //sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]);
          }
          //fprintf(stderr,"\n sam_write works\n");
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->l = 0;
        }
        struct_obj->fqresult_r1->l =0;
        struct_obj->fqresult_r2->l =0;
      }

      memset(qual_r1, 0, sizeof qual_r1); 
      memset(qual_r2, 0, sizeof qual_r2);  
      memset(seq_r1, 0, sizeof seq_r1);
      memset(seq_r1_mod, 0, sizeof seq_r1_mod);

      memset(seq_r2, 0, sizeof seq_r2);
      memset(seq_r2_mod, 0, sizeof seq_r2_mod);
      
      memset(readadapt, 0, sizeof readadapt);
      memset(readadapt2, 0, sizeof readadapt2);
            
      chr_idx = 0;
      //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
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
  /*fprintf(stderr,"\t -> C count %d\n",C_count);
  fprintf(stderr,"\t -> G count %d\n",G_count);
  fprintf(stderr,"\t -> CT count %f\n",CT_count);
  fprintf(stderr,"\t -> CT freq %f\n",CT_count/C_count);
  fprintf(stderr,"\t -> GA count %f\n",GA_count);
  fprintf(stderr,"\t -> GA freq %f\n",GA_count/G_count);*/

  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  //bam_destroy1(struct_obj->list_of_reads[0]);
  for(int j=0; j<struct_obj->m;j++){bam_destroy1(struct_obj->list_of_reads[j]);}

  fprintf(stderr,"\t-> Number of reads generated by thread %d is %d \n",struct_obj->threadno,localread);

  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_flag,
                        const char* Sizefile,int FixedSize,int qualstringoffset,const char* QualProfile1,const char* QualProfile2, int threadwriteno,
                        const char* QualStringFlag,const char* Polynt){
  //creating an array with the arguments to create multiple threads;
  
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];

  // if its defined globally cant i then move nuc2int out of the functions?
  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4; 
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  size_t chr_size_cumm[chr_total+1];
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  size_t genome_size = strlen(genome_data);

  if (genome_data != NULL){
    fprintf(stderr,"\t-> Done creating the large concatenated contig, with size of %lu\n",genome_size);
  
    Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

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
      mode = "ws";
      suffix1 = ".sam";
      alnformatflag++;
    }
    else if(strcasecmp("bam",OutputFormat)==0){
      mode = "wb";
      suffix1 = ".bam";
      alnformatflag++;
    }
    else{fprintf(stderr,"\t-> Fileformat is currently not supported \n");}
    strcat(file1,suffix1);

    /*
    else if(strcasecmp("cram",OutputFormat)==0){
      mode = "wc";
      suffix1 = ".cram"; //Does create a cram file but fix: [E::cram_encode_container] Failed to load reference #0, [E::cram_get_ref] Failed to populate reference for id 0
      alnformatflag++;
    }
    */
    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(alnformatflag == 0){
      fprintf(stderr,"not bam loop \n");
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp1 = bgzf_open(filename1,mode);
      bgzf_mt(bgzf_fp1,mt_cores,bgzf_buf);

      if(strcasecmp("PE",SeqType)==0){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp2 = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp2,mt_cores,bgzf_buf);
      }
      fprintf(stderr,"\t-> Number of cores for bgzf_mt: %d \t %d\n",mt_cores,threadwriteno); 
    }
    else{
      fprintf(stderr,"bam loop \n");
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,seq_ref,chr_total,genome_size);
    }
  
    int number; int* Frag_len; double* Frag_freq;
    if(Sizefile!=NULL){
      Frag_len = new int[LENS];
      Frag_freq = new double[LENS];
      FragArray(number,Frag_len,Frag_freq,Sizefile); //Size_dist_sampling //"Size_dist/Size_freq_modern.txt"
    }
    else{number = -1;}

    // SUPER SAMPLER TEST
    const char *freqfile_r1; //"Qual_profiles/AccFreqL150R1.txt";
    const char *freqfile_r2;
    int outputoffset = qualstringoffset;
    unsigned long readcyclelength;
    ransampl_ws ***QualDist;
    char nt_qual_r1[1024];
    ransampl_ws ***QualDist2;
    char nt_qual_r2[1024];

    double ErrArray_r1[1024];
    double ErrArray_r2[1024];

    const char* bamQflag = QualStringFlag;
    //fprintf(stderr,"\n bamQflag %s\n",QualStringFlag);
    
    if(strcasecmp("true",QualStringFlag)==0){ //|| strcasecmp("bam",OutputFormat)==0
      //fprintf(stderr,"BAMFLAG FOR OUTPUT");
      freqfile_r1 = QualProfile1;
      QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1,readcyclelength);
      if(strcasecmp("PE",SeqType)==0){
        freqfile_r2 = QualProfile2;
        QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2,readcyclelength);
      }
    }
    //fprintf(stderr,"Qualities %c \t %c \t %c \t %c \n",nt_qual_r1[0],nt_qual_r1[1],nt_qual_r1[6],nt_qual_r1[7]);
    //fprintf(stderr,"Err Prob %f \t %f \t %f \t %f \n",ErrArray_r1[0],ErrArray_r1[1],ErrArray_r1[6],ErrArray_r1[7]);

    int maxsize = 20;
    char polynucleotide;
    if (strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}

    //initialzie values that should be used for each thread
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
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;

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

      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (size_t*)malloc(sizeof(size_t) * (struct_for_threads[i].chr_no+1));
      //fprintf(stderr,"struct_for_threads[i]: %p\n",struct_for_threads[i]);
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    
    //fprintf(stderr,"Creating a bunch of threads\n"); 
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Fafq_thread_se_run,&struct_for_threads[i]);
    }
    //fprintf(stderr,"Done Creating a bunch of threads\n");
    

    for (int i = 0; i < nthreads; i++)
    {  
      //fprintf(stderr,"joing threads\n");fflush(stderr);
      pthread_join(mythreads[i],NULL);
      //fprintf(stderr, "\t[ANDET STED] walltime used for join =  %.2f sec\n", (float)(time(NULL) - t3));  
    }
    
    //fprintf(stderr,"Done joining a bunch of threads\n");
    
    if(alnformatflag == 0){
      //fprintf(stderr,"Closing fasta fastq files\n");
      bgzf_close(bgzf_fp1);
      if(strcasecmp("PE",SeqType)==0){bgzf_close(bgzf_fp2);}
    }
    else{
      sam_hdr_destroy(SAMHeader);
      sam_close(SAMout);
    } 
    
    //fprintf(stderr,"Before memory cleaning\n");

    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      free(struct_for_threads[i].list_of_reads);
      delete struct_for_threads[i].fqresult_r1;

      free(struct_for_threads[i].fqresult_r2 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      delete struct_for_threads[i].fqresult_r2;      
    }
    

    if(strcasecmp("true",QualStringFlag)==0){ //if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0 || bamQflag == "TRUE"){ //strcasecmp("bam",OutputFormat)==0
      //fprintf(stderr,"INSIDE TRUE QUALSTRING IF FINAL\n");
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
    
    free(genome_data);
    fflush(stderr);
  }
  return NULL;
}
//./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 1000000 -s 1 -t1 2 -lf Size_dist/Size_dist_sampling.txt -b 0.024,0.36,0.68,0.0097 -seq SE -f bam -o chr22