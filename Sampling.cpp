#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>

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
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "ThreadGeneration.h"
#include "Sampling.h"
#include "sample_qscores.h"
#include "fasta_sampler.h"

#define LENS 4096
#define MAXBINS 100
extern int refToInt[256];

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void ReversComplement(char seq[]){
  // generates the reverse complementary sequence from an input sequence

  char NtComp[5] = {'T', 'G', 'C', 'A','N'};
  char seq_intermediate[1024] = {0};
  strcpy(seq_intermediate,seq);
  //fprintf(stderr,"SEQUENCE \t\t%s\n",seq_intermediate);
  int seqlen = strlen(seq);
  //Complementing sequence
  for(int i=0;i<seqlen;i++){
    seq_intermediate[i] = NtComp[refToInt[(unsigned char) seq_intermediate[i]]]; //warning: array subscript has type 'char' [-Wchar-subscripts]
  }
  //fprintf(stderr,"COMP SEQUENCE \t\t%s\n",seq_intermediate);

  //reverse complement
  for(int i=seqlen-1;i>-1;i--){
    seq[seqlen-i-1] = seq_intermediate[i];
  }
  //just to ensure no issues arise in case of not clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}


pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  mrand_t *drand_alloc = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt_adapt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_briggs = mrand_alloc(struct_obj->rng_type,loc_seed);

  char *seq;//actual sequence, this is unallocated
  int fraglength;

  // sequence reads, original, modified, with adapters, pe
  char READ_ID[1024];
  char seq_r1[1024] = {0};
  char seq_r2[1024] = {0};
  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};
 
  
  size_t reads = struct_obj -> reads;
  size_t BufferLength = struct_obj -> BufferLength;

  int ErrProbTypeOffset = 0;
  if (struct_obj->OutputFormat==fqT || struct_obj->OutputFormat==fqgzT){ErrProbTypeOffset=33;}
  kstring_t *fqs[2];
  for(int i=0;i<2;i++){
    fqs[i] =(kstring_t*) calloc(1,sizeof(kstring_t));
    fqs[i]->s = NULL;
    fqs[i]->l = fqs[i]->m = 0;
  }
  
  size_t localread = 0;
  int iter = 0;
  size_t current_reads_atom = 0;

  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int posB,posE;//this is the first and last position of our fragment

  extern int SIG_COND;
<<<<<<< HEAD
  while (current_reads_atom < reads &&SIG_COND) {
    double rand_val = mrand_pop(drand_alloc);
    rand_start = rand_val * (genome_len-300)+1; //genome_len-100000;
    double rand_val_len = mrand_pop(drand_alloc);
    // Fragment length creation
    int fraglength;
    if (struct_obj->No_Len_Val != -1){
      // random start and length are not dependent on the same rand val
      int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val_len);
      fraglength =  struct_obj->FragLen[lengthbin];
    }
    else{
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
    int rand_id = (rand_val_id * fraglength-1); //100
    
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}
    if (fraglength > readsizelimit){strncpy(seq_r1,struct_obj->genome+rand_start-1,readsizelimit);}   // case 1
    else {strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    
    if(strcasecmp("PE",struct_obj->SeqType)==0){
      if (fraglength > readsizelimit){
=======
>>>>>>> development

  std::default_random_engine RndGen(loc_seed);
  
  sim_fragment *sf;
  if (struct_obj->LengthType==0)
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->FixedSize,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  else
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->distparam1,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  
  size_t moduloread = reads/10;
  while (current_reads_atom < reads && SIG_COND){
    //lets start by resetting out datastructures to NULL, nill, nothing.
    qual_r1[0] = qual_r2[0] = seq_r1[0] = seq_r2[0] = '\0';

<<<<<<< HEAD
    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N'); //first encounter with N
    pch2 = strrchr(seq_r1,'N'); //last encounter with N
=======
    //sample fragmentlength
    int fraglength = getFragmentLength(sf); //fraglength = abs(mrand_pop_long(drand_alloc)) % 1000;
    
    // Selecting genomic start position across the generated contiguous contigs for which to extract 
    int chr_idx = -1;

    //maxbases is the number of nucleotides we will work with.
    //set this to minimum of bases sequenced or fragment length.
    //if output is fasta then legnth is simply the entire fragment
    int maxbases = std::min(fraglength,struct_obj->maxreadlength);
    if (struct_obj->OutputFormat==faT ||struct_obj->OutputFormat==fagzT)
      maxbases = fraglength;

    //Generating random ID unique for each read output
    double rand_val_id = mrand_pop(drand_alloc);
    //fprintf(stderr,"RAndom values for ID %lf \t %lf\n",rand_val_id,mrand_pop(drand_alloc));
    int rand_id = (rand_val_id * fraglength-1)+(rand_val_id*current_reads_atom); //100
    //why above? just to get random unique?

    //get shallow copy of chromosome, offset into, is defined by posB, and posE
    char *seq = sample(struct_obj->reffasta,drand_alloc,&chr,chr_idx,posB,posE,fraglength);

    //fprintf(stderr,"chr %s \t idx %d \n",chr,chr_idx);
    //now copy the actual sequence into seq_r1 and seq_r2 if PE 
    strncpy(seq_r1,seq+posB,maxbases);
    if(PE==struct_obj->SeqType)
      strncpy(seq_r2,seq+posE-maxbases,maxbases);
    
    //fprintf(stderr,"CHECKING THE LENGTH ISSUES %d \t %d \t %d \t %d \n",fraglength,struct_obj->maxreadlength,maxbases,strlen(seq_r1));
>>>>>>> development

    /*
      ||------R1------>,,,,,,,,,,|-------R2----->||
    */
    
    //fprintf(stderr,"Current read %zu and max bases %d and fraglent %d  and coordinates %d \t %d \t seq length %d\n",current_reads_atom,maxbases,fraglength,posB,posE,strlen(seq_r1));
    //Selecting strand, 0-> forward strand (+) 5'->3', 1 -> reverse strand (-) 3'->5'
    //rename to strand to strandR1
    int strandR1 = mrand_pop(drand_alloc)>0.5?0:1; //int strand = (int) (rand_start%2);
    
    //Remove reads which are all N? In this code we remove both pairs if both reads are all N
    int skipread = 1;
    for(int i=0;skipread&&i<(int)strlen(seq_r1);i++)
      if(seq_r1[i]!='N')
	      skipread = 0;

    for(int i=0; skipread && i<(int)strlen(seq_r2);i++)
      if(seq_r2[i]!='N')
	    skipread = 0;
    
    if(skipread==1)
      continue;
    
    // generating sam output information
    int seqlen = strlen(seq_r1);

    if(strlen(seq_r1) < 20)
      continue;

    int SamFlags[2] = {-1,-1}; //flag[0] is for read1, flag[1] is for read2
    
    //now everything is the same strand as reference, which we call plus/+
    //lets flip to 5 to 3
    if (SE==struct_obj->SeqType){
      if (strandR1 == 0)
	      SamFlags[0] = 0;
      else if (strandR1 == 1){
        SamFlags[0] = 16;
        ReversComplement(seq_r1);
      }
    }
    if (PE==struct_obj->SeqType){
      if (strandR1 == 0){
        SamFlags[0] = 97;
        SamFlags[1] = 145;
        ReversComplement(seq_r2);
      }
      else if (strandR1 == 1){
        SamFlags[0] = 81;
        SamFlags[1] = 161;
        ReversComplement(seq_r1);
      }
    }
<<<<<<< HEAD
    else{
      // then all the same start pos would have the same strand no matter the chromosome?
      int strand = (int) (rand_start%2);//(int) rand_r(&loc_seed)%2;//1;//rand() % 2;
      
      if (struct_obj->SAMout){
        // in bam files the reads are all aligning to the forward strand, but the flags identify read property
        if (strcasecmp("SE",struct_obj->SeqType)==0){
          if (strand == 0){flag = 0;} // se read forward strand
          else{flag = 16;
            // Simulate reads from the negative strand but change the orientation to 5' to 3' fragment
            DNA_complement(seq_r1);
            reverseChar(seq_r1,strlen(seq_r1));
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
=======
    //so now everything is on 5->3 and some of them will be reverse complement to referene

    //Nucleotide alteration models only on the sequence itself which holds for fa,fq,sam
    if(struct_obj->DoBriggs){
      SimBriggsModel(seq_r1,fraglength,struct_obj->BriggsParam[0],
		     struct_obj->BriggsParam[1],
		     struct_obj->BriggsParam[2], 
		     struct_obj->BriggsParam[3],drand_alloc_briggs);
      if (PE==struct_obj->SeqType){
        SimBriggsModel(seq_r2, fraglength,struct_obj->BriggsParam[0], 
		       struct_obj->BriggsParam[1], 
		       struct_obj->BriggsParam[2], 
		       struct_obj->BriggsParam[3],drand_alloc_briggs);
>>>>>>> development
      }
    } 

<<<<<<< HEAD
        if(strcasecmp("true",struct_obj->SubFlag)==0){
          MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          if (strcasecmp("PE",struct_obj->SeqType)==0){
            MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
          }
        }
      }
=======
    //should this be based on forward + strand or on the 5->3 sequences?
    if(struct_obj->doMisMatchErr){
      //this function below modifies sequence to include substitutions from mismatch incorperation file
      //this option is currently -mf, NB change name
      MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
      if (PE==struct_obj->SeqType)
	      MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
    }
>>>>>>> development
      
    
    sprintf(READ_ID,"T%d_RID%d_S%d_%s:%d-%d_length:%d", struct_obj->threadno, rand_id,strandR1,chr,posB+1,posE,fraglength);
    
    int nsofts[2] = {0,0};//this will contain the softclip information to be used by sam/bam/cram out
    //below will contain the number of bases for R1 and R2 that should align to reference before adding adapters and polytail
    int naligned[2] = {(int)strlen(seq_r1),-1};
    if(PE==struct_obj->SeqType)
      naligned[1] = strlen(seq_r2);
    
    
    //add adapters
    if(struct_obj->AddAdapt){
      // Because i have reverse complemented the correct sequences and adapters depending on the strand origin (or flags), i know all adapters will be in 3' end
      nsofts[0] = std::min(struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
      //fprintf(stderr,"The minimum values are %d \t %d \t %d \n",maxbases,struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
      strncpy(seq_r1+strlen(seq_r1),struct_obj->Adapter_1,nsofts[0]);
      if(PE==struct_obj->SeqType){
        nsofts[1] = std::min(struct_obj->maxreadlength-strlen(seq_r2),strlen(struct_obj->Adapter_2));
        strncpy(seq_r2+strlen(seq_r2),struct_obj->Adapter_2,nsofts[1]);
      }
    }

    //add polytail
    if (struct_obj->PolyNt != 'F') {
      int nitems = struct_obj->maxreadlength-strlen(seq_r1);
      memset(seq_r1+strlen(seq_r1),struct_obj->PolyNt,nitems);
      nsofts[0] += nitems;
      if(PE==struct_obj->SeqType){
        nitems = struct_obj->maxreadlength-strlen(seq_r2);
        memset(seq_r2+strlen(seq_r2),struct_obj->PolyNt,nitems);
        nsofts[1] += nitems;
      }
    }

    if(struct_obj->SAMout){
      //sanity check
      //fprintf(stderr,"SANITY CHECK seq_R1 %d \t %d \t %d \n seq_R2  %d \t %d \t %d \n",strlen(seq_r1),naligned[0],nsofts[0],strlen(seq_r2),naligned[1],nsofts[1]);
      if(strlen(seq_r1)!=naligned[0]+nsofts[0]){
        fprintf(stderr,"Number of aligned bases + number of adap + poly does not match\n");
        exit(0);
      }
      //below only runs for PE that is when nalign[1] is not -1
      if(naligned[1]!=-1 && strlen(seq_r2)!=naligned[1]+nsofts[1]){
        fprintf(stderr,"Number of aligned bases + number of adap + poly does not match\n");
        exit(0);
      }
    }
    //now seq_r1 and seq_r2 is completely populated let is build qualscore if
    //saving both fasta and adapter to fasta format
    if (struct_obj->OutputFormat==faT ||struct_obj->OutputFormat==fagzT){
      ksprintf(fqs[0],">%s R1\n%s\n",READ_ID,seq_r1);//make this into read
      if (PE==struct_obj->SeqType)
	    ksprintf(fqs[1],">%s R2\n%s\n",READ_ID,seq_r2);
    } 
    else{
      // Fastq and Sam needs quality scores
      //if(strandR1==0){fprintf(stderr,"----------\nSEQUENCE \n%s\n",seq_r1);}//int ntcharoffset
      sample_qscores(seq_r1,qual_r1,strlen(seq_r1),struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt_adapt,struct_obj->DoSeqErr,ErrProbTypeOffset);
      //if(strandR1==0){fprintf(stderr,"%s\n",seq_r1);}
      if (PE==struct_obj->SeqType)
      	sample_qscores(seq_r2,qual_r2,strlen(seq_r2),struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt_adapt,struct_obj->DoSeqErr,ErrProbTypeOffset);
      
      //write fq if requested
      if (struct_obj->OutputFormat==fqT || struct_obj->OutputFormat==fqgzT){
        ksprintf(fqs[0],"@%s R1\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
        if (PE==struct_obj->SeqType)
          ksprintf(fqs[1],"@%s R2\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);
      }

      //now only sam family needs to be done, lets revcomplement the bases, and reverse the quality scores so everything is back to forward/+ strand
      if(struct_obj->SAMout){
        uint32_t AlignCigar[2][10];//cigs[0] is read1 cigs[1] is read2
        size_t n_cigar[2] = {1,1};

      if (struct_obj->NoAlign){
        AlignCigar[0][0] = bam_cigar_gen(naligned[0], BAM_CMATCH);
        if(nsofts[0]>0){
          AlignCigar[0][1] = bam_cigar_gen(nsofts[0], BAM_CSOFT_CLIP);
          n_cigar[0] = 2;
        }
<<<<<<< HEAD
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          for(long unsigned int p = 0;p<strlen(readadapt);p++){
            //double rand_val = mrand_pop(drand_alloc);
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
            sprintf(readadapt_rc_sam1, "%s", read_rc_sam1);
            strcat(readadapt_rc_sam1, readadapt_err1);
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
            
            sprintf(readadapt_rc_sam2, "%s", read_rc_sam2);
            strcat(readadapt_rc_sam2, readadapt_err2);
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
=======
        if(strandR1==1){
          ReversComplement(seq_r1);
          reverseChar(qual_r1,strlen(seq_r1));
          if(n_cigar[0]>1){
            //swap softclip and match
            uint32_t tmp= AlignCigar[0][0];
            AlignCigar[0][0] = AlignCigar[0][1];
            AlignCigar[0][1] = tmp;
          }
        }
        if (PE==struct_obj->SeqType){
          AlignCigar[1][0] = bam_cigar_gen(naligned[1], BAM_CMATCH);
          if(nsofts[1]>0){
            AlignCigar[1][1] = bam_cigar_gen(nsofts[1], BAM_CSOFT_CLIP);
            n_cigar[1] = 2;
          }
          if(strandR1==0){
            ReversComplement(seq_r2);
            reverseChar(qual_r2,strlen(seq_r2));
            if(n_cigar[1]>1){
              //swap softclip and match
              uint32_t tmp= AlignCigar[1][0];
              AlignCigar[1][0] = AlignCigar[1][1];
              AlignCigar[1][1] = tmp;
>>>>>>> development
            } 
          }
        }
      }
      else{
<<<<<<< HEAD
        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s_R1\n%s\n",READ_ID,seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s_R2\n%s\n",READ_ID,seq_r2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
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
          if (flag == 16 || flag == 81){DNA_complement(seq_r1);reverseChar(seq_r1,strlen(seq_r1));}
          else if (flag == 97){DNA_complement(seq_r2);reverseChar(seq_r2,strlen(seq_r2));}  
          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
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
=======
        //this is unaligned part
        AlignCigar[0][0] = bam_cigar_gen(strlen(seq_r1), BAM_CSOFT_CLIP);
        AlignCigar[1][0] = bam_cigar_gen(strlen(seq_r2), BAM_CSOFT_CLIP);
>>>>>>> development
      }
      //now reads, cigards and quals are correct 

      //generating id, position and the remaining sam field information
      size_t l_aux = 2; uint8_t mapq = 60;
      hts_pos_t min_beg, max_end, insert; //max_end, insert;
      min_beg = posB;
      max_end = posE;
      hts_pos_t min_beg_mate, max_end_mate;
      insert = max_end - min_beg;
      /*if (PE==struct_obj->SeqType)
        insert = max_end - min_beg + 1;
      else
        insert = max_end - min_beg + 1;
        //insert = max_end = min_beg = -1;*/

<<<<<<< HEAD
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
=======
      const char* suffR1 = " R1";
      const char* suffR2 = " R2";
>>>>>>> development

      char READ_ID2[1024];
      strcpy(READ_ID2,READ_ID);
      
      strcat(READ_ID,suffR1);
      //fprintf(stderr,"CHR IDX %d\n",chr_idx);
      //fprintf(stderr,"BEGIN %d AND END %d\n",min_beg,max_end);
      int chr_idx_mate = -1;
      int chr_max_end_mate = -1;
      int insert_mate = 0;
      if(PE==struct_obj->SeqType){
        strcat(READ_ID2,suffR2);
        chr_idx_mate = chr_idx;
        chr_max_end_mate = max_end;
        insert_mate = insert;
        //fprintf(stderr,"CHR IDX %d\n",chr_idx_mate);
      }
      if (struct_obj->NoAlign == 0){
        mapq = 255;
        SamFlags[0] = SamFlags[1] = 4;
        chr_idx = -1;
        min_beg = max_end = -1;
      }

<<<<<<< HEAD
        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) != assert_int);
          }
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->l = 0;
=======
      //we have set the parameters accordingly above for no align and PE
      //fprintf(stderr,"chr idx %d \t chr_max %d \t chr_insert %d\n",chr_idx_mate,chr_max_end_mate,insert_mate);
      bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID),READ_ID,SamFlags[0],chr_idx,min_beg,mapq,
            n_cigar[0],AlignCigar[0],chr_idx_mate,chr_max_end_mate-1,insert_mate,strlen(seq_r1),seq_r1,qual_r1,l_aux);
      //exit(0);
      //write PE also
      if (PE==struct_obj->SeqType){
        bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID2),READ_ID2,SamFlags[1],chr_idx,max_end-1,mapq,
          n_cigar[1],AlignCigar[1],chr_idx,min_beg,0-insert_mate,strlen(seq_r2),seq_r2,qual_r2,l_aux);
      }
    
      if (struct_obj->LengthData < struct_obj->MaximumLength){   
        pthread_mutex_lock(&write_mutex);
        for (int k = 0; k < struct_obj->LengthData; k++){
          assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) >=0 );
>>>>>>> development
        }
        pthread_mutex_unlock(&write_mutex);
        struct_obj->LengthData = 0;
      }
      fqs[0]->l = fqs[1]->l = 0;
      }
      
<<<<<<< HEAD
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
=======
>>>>>>> development
    }
    
    if (struct_obj->bgzf_fp[0]){
      if (fqs[0]->l > BufferLength){
        pthread_mutex_lock(&write_mutex);
        assert(bgzf_write(struct_obj->bgzf_fp[0],fqs[0]->s,fqs[0]->l)!=0);
        if (PE==struct_obj->SeqType){
          assert(bgzf_write(struct_obj->bgzf_fp[1],fqs[1]->s,fqs[1]->l)!=0);
        }
        pthread_mutex_unlock(&write_mutex);
	fqs[0]->l = fqs[1]->l = 0;
      }
    }

    memset(qual_r1, 0, sizeof qual_r1); 
    memset(qual_r2, 0, sizeof qual_r2);  
    memset(seq_r1, 0, sizeof seq_r1);
    memset(seq_r2, 0, sizeof seq_r2);

    chr_idx = -1;
    iter++;
    localread++;
    current_reads_atom++;
    //printing out every tenth of the runtime
    if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
      fprintf(stderr,"\t-> Thread %d produced %zu reads with a current total of %zu\n",struct_obj->threadno,moduloread,current_reads_atom);
  }
  if (struct_obj->bgzf_fp[0]){
    if (fqs[0]->l > 0){
      pthread_mutex_lock(&write_mutex);
      assert(bgzf_write(struct_obj->bgzf_fp[0],fqs[0]->s,fqs[0]->l)!=0);
      if (PE==struct_obj->SeqType){assert(bgzf_write(struct_obj->bgzf_fp[1],fqs[1]->s,fqs[1]->l)!=0);}
      pthread_mutex_unlock(&write_mutex);
	fqs[0]->l = fqs[1]->l = 0;
    } 
  }

  for(int j=0; j<struct_obj->MaximumLength;j++){bam_destroy1(struct_obj->list_of_reads[j]);}
  
  free(sf->rand_alloc);
  delete sf;

  if(fqs[0]->s)
    free(fqs[0]->s);
  if(fqs[1]->s)
    free(fqs[1]->s);

  free(fqs[0]);
  free(fqs[1]);
  
  free(drand_alloc);
  free(drand_alloc_nt);
  free(drand_alloc_nt_adapt);
  free(drand_alloc_briggs);
  
  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);

<<<<<<< HEAD
  pthread_exit(NULL);
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, size_t reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_flag,
                        const char* Sizefile,int FixedSize,int SizeDistType, int val1, int val2,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2, int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,const char* ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* MisMatchFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat,char* Variant_flag,const char *VarType,
                        char CommandArray[1024],const char* version,const char* HeaderIndiv,const char* NoAlign,size_t BufferLength){
  //creating an array with the arguments to create multiple threads;

  int nthreads=thread_no;
  pthread_t mythreads[nthreads];

  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4; 

  int chr_total = 0;
  char *genome_data;
  if (Specific_Chr[0] != NULL){
    while (Specific_Chr[chr_total]){
      chr_total++;
      }
  }
  else{
    chr_total = faidx_nseq(seq_ref);
  }

  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_idx_arr[chr_total];
  size_t chr_size_cumm[chr_total+1];
  
  if (chr_total < faidx_nseq(seq_ref)){
    for (int j = 0; j < faidx_nseq(seq_ref); j++){
      for (int i = 0; i < chr_total; i++){
        if(strcasecmp(faidx_iseq(seq_ref, j),Specific_Chr[i])==0){
          chr_idx_arr[i] = j;
        }
      } 
    }
  }
  else
  {
    for (int j = 0; j < faidx_nseq(seq_ref); j++){;chr_idx_arr[j] = j;}
  }
  
  if(VCFformat != NULL && strcasecmp(Variant_flag,"bcf")==0){
    //const char* HeaderIndiv = "HG00097";
    genome_data = full_vcf_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm,VCFformat,VarType,HeaderIndiv);
  }
  else{
    if (chr_total == faidx_nseq(seq_ref)){
      genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
    }  
    else{
      genome_data = partial_genome_create(seq_ref,chr_total-1,chr_sizes,Specific_Chr,chr_size_cumm);
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
    htsThreadPool p = {NULL, 0};

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
      else{suffix1 = "_R1.fa";suffix2 = "_R2.fa";}
    }
    else if(strcasecmp("fa.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa.gz";}
      else{suffix1 = "_R1.fa.gz";suffix2 = "_R2.fa.gz";}
    }
    else if(strcasecmp("fq",OutputFormat)==0){
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq";}
      else{suffix1 = "_R1.fq";suffix2 = "_R2.fq";}
    }
    else if(strcasecmp("fq.gz",OutputFormat)==0){
      mode = "w";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq.gz";}
      else{suffix1 = "_R1.fq.gz";suffix2 = "_R2.fq.gz";}
    }
    else if(strcasecmp("sam",OutputFormat)==0){
      mode = "ws";
      suffix1 = ".sam";
      alnformatflag++;
    }
    else if(strcasecmp("bam",OutputFormat)==0){
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
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp1 = bgzf_open(filename1,mode); //w
      bgzf_mt(bgzf_fp1,mt_cores,bgzf_buf); //
      
      if(strcasecmp("PE",SeqType)==0){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp2 = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp2,mt_cores,bgzf_buf);
      }
    }
    else{
      char *ref =(char*) malloc(10 + strlen(FastaFileName) + 1);
      sprintf(ref, "reference=%s", FastaFileName);
      
      hts_opt_add((hts_opt **)&fmt_hts->specific,ref);
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();

      if(threadwriteno>0){
        if (!(p.pool = hts_tpool_init(threadwriteno))) {
          fprintf(stderr, "Error creating thread pool\n");
          exit(0);
        }
        hts_set_opt(SAMout, HTS_OPT_THREAD_POOL, &p);
      }
      // generate header
      //hts_set_threads(SAMout, 4);
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,seq_ref,chr_total,chr_idx_arr,genome_size,CommandArray,version);
      free(ref);
      hts_opt_free((hts_opt *)fmt_hts->specific);
    }

    int number; int* Frag_len; double* Frag_freq;
    if(FixedSize==-1 && SizeDistType==-1){
      Frag_len = new int[LENS];
      Frag_freq = new double[LENS];
      FragArray(number,Frag_len,Frag_freq,Sizefile); //Size_dist_sampling //"Size_dist/Size_freq_modern.txt"
    }
    else if(FixedSize==-1 && SizeDistType!=-1){
      Frag_len = new int[LENS];
      Frag_freq = new double[LENS];
      FragDistArray(number,Frag_len,Frag_freq,SizeDistType,seed,val1, val2);
    }
    else if(FixedSize!=-1){number = -1;}
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
    if(strcasecmp("true",QualStringFlag)==0){ //|| strcasecmp("bam",OutputFormat)==0
      freqfile_r1 = QualProfile1;
      QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1,readcyclelength);
      if(strcasecmp("PE",SeqType)==0){
        freqfile_r2 = QualProfile2;
        QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2,readcyclelength);
      }
    }
    
    int maxsize = 20;
    char polynucleotide;
    if (Polynt != NULL && strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}

    double* MisMatchFreqArray;
    int mismatchcyclelength = 0;
    if (SubProfile != NULL){
      MisMatchFreqArray = new double[LENS];
      MisMatchFreqArray = MisMatchFileArray(MisMatchFreqArray,SubProfile,mismatchcyclelength);
    }

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

      struct_for_threads[i].MisMatch = MisMatchFreqArray;
      struct_for_threads[i].SubFlag = MisMatchFlag;
      struct_for_threads[i].MisLength = (int) mismatchcyclelength;
      struct_for_threads[i].readcycle = (int) readcyclelength;
      struct_for_threads[i].reads = reads;
      struct_for_threads[i].BufferLength = BufferLength;

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
      struct_for_threads[i].NoAlign = (char) NoAlign[0];
      struct_for_threads[i].Variant_flag = Variant_flag;


      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (size_t*)malloc(sizeof(size_t) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (char**)malloc(sizeof(char*) * struct_for_threads[i].chr_no+1);
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
    if(FixedSize == -1){
      delete[] Frag_freq;
      delete[] Frag_len;
    }
    
    if(SubProfile != NULL){delete[] MisMatchFreqArray;}
    
    free(genome_data);
    fflush(stderr);
  }
  return NULL;
=======
  if(struct_obj->totalThreads>1)
    pthread_exit(NULL);
  
  return 0;
>>>>>>> development
}
