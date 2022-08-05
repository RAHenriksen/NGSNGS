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
#include "sample_qscores.h"
#include "fasta_sampler.h"

#define LENS 4096
#define MAXBINS 100
unsigned char nuc2int[255];

pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  mrand_t *drand_alloc = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt_adapt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_briggs = mrand_alloc(struct_obj->rng_type,loc_seed);

  fprintf(stderr,"INSIDE SAMPLING THREADS\n");
  fprintf(stderr,"-----------------------\nFASTA SAMPLER \n---------------------\n");

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

  
  size_t localread = 0;
  int iter = 0;
  size_t current_reads_atom = 0;

  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int posB,posE;//this is the first and last position of our fragment

  extern int SIG_COND;

  std::default_random_engine RndGen(loc_seed);
  
  sim_fragment *sf;
  if (struct_obj->LengthType==0)
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->FixedSize,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  else
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->distparam1,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  
  size_t moduloread = reads/10;
  while (current_reads_atom < reads && SIG_COND){
    //lets start by resetting out datastructures to NULL, nill, nothing.
    memset(qual_r1, 0, sizeof qual_r1); 
    memset(qual_r2, 0, sizeof qual_r2);  

    memset(seq_r1, 0, sizeof seq_r1);
    memset(seq_r2, 0, sizeof seq_r2);
      
    //printing out every tenth of the runtime
    if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
      fprintf(stderr,"\t-> Thread %d roduced %zu reads with a current total of %zu\n",struct_obj->threadno,moduloread,current_reads_atom);

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
    int rand_id = (rand_val_id * fraglength-1); //100
    //why above? just to get random unique?

    //get shallow copy of chromosome, offset into, is defined by posB, and posE
    char *seq = sample(struct_obj->reffasta,drand_alloc,&chr,chr_idx,posB,posE,fraglength);

    //now copy the actual sequence into seq_r1 and seq_r2 if PE 
    strncpy(seq_r1,seq+posB,maxbases);
    if(PE==struct_obj->SeqType)
      strncpy(seq_r2,seq+posE-maxbases,maxbases);
    
    /*
      ||------R1------>,,,,,,,,,,|-------R2----->||
    */
    
    //Selecting strand, 0-> forward strand (+) 5'->3', 1 -> reverse strand (-) 3'->5'
    //rename to strand to strandR1
    int strandR1 = mrand_pop(drand_alloc)>0.5?0:1; //int strand = (int) (rand_start%2);
    
    //Remove reads which are all N? In this code we remove both pairs if both reads are all N
    int skipread = 1;
    for(int i=0;skipread&&i<(int)strlen(seq_r1);i++)
      if(seq_r1[i]!='N')
	      skipread = 0;

    for(int i=0;seq_r2 && skipread && i<(int)strlen(seq_r2);i++)
      if(seq_r2[i]!='N')
	    skipread = 0;
    
    if(skipread==1)
      continue;
    
    // generating sam output information
    int seqlen = strlen(seq_r1);
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
    //so now everything is on 5->3 and some of them will be reverse complement to referene

    //Nucleotide alteration models only on the sequence itself which holds for fa,fq,sam
    if(struct_obj->DoBriggs){
      SimBriggsModel(seq_r1,fraglength,struct_obj->BriggsParam[0],
		     struct_obj->BriggsParam[1],
		     struct_obj->BriggsParam[2], 
		     struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
      if (PE==struct_obj->SeqType){
        SimBriggsModel(seq_r2, fraglength,struct_obj->BriggsParam[0], 
		       struct_obj->BriggsParam[1], 
		       struct_obj->BriggsParam[2], 
		       struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
      }
    }

    //should this be based on forward + strand or on the 5->3 sequences?
    if(struct_obj->doMisMatchErr){
      //this function below modifies sequence to include substitutions from mismatch incorperation file
      //this option is currently -mf, NB change name
      MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
      if (PE==struct_obj->SeqType)
	      MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
    }
      
    
    sprintf(READ_ID,"T%d_RID%d_S%d_%s:%zu-%zu_length:%d", struct_obj->threadno, rand_id,strandR1,chr,posB,posE,fraglength);
    
    int nsofts[2] = {0,0};//this will contain the softclip information to be used by sam/bam/cram out
    //below will contain the number of bases for R1 and R2 that should align to reference before adding adapters and polytail
    int naligned[2] = {(int)strlen(seq_r1),-1};
    if(PE==struct_obj->SeqType)
      naligned[1] = strlen(seq_r2);
    
    
    //add adapters
    if(struct_obj->AddAdapt){
      //fprintf(stderr,"INSIDE ADD ADAPT!\n%s\n",seq_r1);
      // Because i have reverse complemented the correct sequences and adapters depending on the strand origin (or flags), i know all adapters will be in 3' end
      nsofts[0] = std::min(struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
      //fprintf(stderr,"The minimum values are %d \t %d \t %d \n",maxbases,struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
      strncpy(seq_r1+strlen(seq_r1),struct_obj->Adapter_1,nsofts[0]);
      //fprintf(stderr,"SEQUENCE R1\n%s\n",seq_r1);
      if(PE==struct_obj->SeqType){
        nsofts[1] = std::min(struct_obj->maxreadlength-strlen(seq_r2),strlen(struct_obj->Adapter_2));
        //fprintf(stderr,"The minimum values are %d \t %d \t %d \n",struct_obj->maxreadlength,struct_obj->maxreadlength-strlen(seq_r2),strlen(struct_obj->Adapter_2));
        //fprintf(stderr,"INSIDE ADD ADAPT!\n%s\n",seq_r2);
        strncpy(seq_r2+strlen(seq_r2),struct_obj->Adapter_2,nsofts[1]);
        //fprintf(stderr,"SEQUENCE R1\n%s\n",seq_r2);
      }
      //do we need to keep track of both or is there some kind of symmetry?
    }

    //add polytail
    if (struct_obj->PolyNt != 'F') {
      int nitems = struct_obj->maxreadlength-strlen(seq_r1);
      //fprintf(stderr,"The minimum values are %d \t %d \t %d \n",struct_obj->maxreadlength,strlen(seq_r1),nitems);
      //fprintf(stderr,"INSIDE POLY ADAPT!\n%s\n",seq_r1);
      memset(seq_r1+strlen(seq_r1),struct_obj->PolyNt,nitems);
      //fprintf(stderr,"INSIDE POLY ADAPT!\n%s\n",seq_r1);
      //fprintf(stderr,"SOFT CLIPPED NUMEBR\n%d\n",nsofts[0]);
      nsofts[0] += nitems;
      //fprintf(stderr,"SOFT CLIPPED NUMEBR\n%d\n",nsofts[0]);
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
      ksprintf(struct_obj->fqresult_r1,">%s R1\n%s\n",READ_ID,seq_r1);//make this into read
      if (PE==struct_obj->SeqType)
	      ksprintf(struct_obj->fqresult_r2,">%s R2\n%s\n",READ_ID,seq_r2);
    } 
    else{
      // Fastq and Sam needs quality scores
      sample_qscores(seq_r1,qual_r1,strlen(seq_r1),struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt_adapt,struct_obj->DoSeqErr);
      if (PE==struct_obj->SeqType)
      	sample_qscores(seq_r2,qual_r2,strlen(seq_r2),struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt_adapt,struct_obj->DoSeqErr);
      
      //write fq if requested
      if (struct_obj->OutputFormat==fqT || struct_obj->OutputFormat==fqgzT){
        ksprintf(struct_obj->fqresult_r1,"@%s R1\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
        if (PE==struct_obj->SeqType)
          ksprintf(struct_obj->fqresult_r2,"@%s R2\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);
      }

      //now only sam family needs to be done, lets revcomplement the bases, and reverse the quality scores so everything is back to forward/+ strand
      if(struct_obj->SAMout){
        uint32_t AlignCigar[2][10];//cigs[0] is read1 cigs[1] is read2
        size_t n_cigar[2] = {1,1};

      if (struct_obj->NoAlign != 'T'){
        AlignCigar[0][0] = bam_cigar_gen(naligned[0], BAM_CMATCH);
        if(nsofts[0]>0){
          AlignCigar[0][1] = bam_cigar_gen(nsofts[0], BAM_CSOFT_CLIP);
          n_cigar[0] = 2;
        }
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
            } 
          }
        }
      }
      else{
        //this is unaligned part
        AlignCigar[0][0] = bam_cigar_gen(strlen(seq_r1), BAM_CSOFT_CLIP);
        AlignCigar[1][0] = bam_cigar_gen(strlen(seq_r2), BAM_CSOFT_CLIP);
      }
      //now reads, cigards and quals are correct 

      //generating id, position and the remaining sam field information
      size_t l_aux = 2; uint8_t mapq = 60;
      hts_pos_t min_beg, max_end, insert; //max_end, insert;
      min_beg = posB;
      max_end = posE;
      hts_pos_t min_beg_mate, max_end_mate;
      insert = max_end - min_beg + 1;
      /*if (PE==struct_obj->SeqType)
        insert = max_end - min_beg + 1;
      else
        insert = max_end - min_beg + 1;
        //insert = max_end = min_beg = -1;*/

      const char* suffR1 = " R1";
      const char* suffR2 = " R2";

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
      if (struct_obj->NoAlign == 'T'){
        mapq = 255;
        SamFlags[0] = SamFlags[1] = 4;
        chr_idx = -1;
        min_beg = max_end = -1;
      }

      //we have set the parameters accordingly above for no align and PE
      //fprintf(stderr,"chr idx %d \t chr_max %d \t chr_insert %d\n",chr_idx_mate,chr_max_end_mate,insert_mate);
      bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID),READ_ID,SamFlags[0],chr_idx,min_beg,mapq,
            n_cigar[0],AlignCigar[0],chr_idx_mate,chr_max_end_mate,insert_mate,strlen(seq_r1),seq_r1,qual_r1,l_aux);
      //exit(0);
      //write PE also
      if (PE==struct_obj->SeqType){
        bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID2),READ_ID2,SamFlags[1],chr_idx,max_end,mapq,
          n_cigar[1],AlignCigar[1],chr_idx,min_beg,0-insert_mate,strlen(seq_r2),seq_r2,qual_r2,l_aux);
      }
    
      if (struct_obj->LengthData < struct_obj->MaximumLength){   
        pthread_mutex_lock(&write_mutex);
        for (int k = 0; k < struct_obj->LengthData; k++){
          assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) >=0 );
        }
        pthread_mutex_unlock(&write_mutex);
        struct_obj->LengthData = 0;
      }
      struct_obj->fqresult_r1->l =0;
      struct_obj->fqresult_r2->l =0;	
      }
      
    }
    
    if (struct_obj->bgzf_fp[0]){
      if (struct_obj->fqresult_r1->l > BufferLength){
        pthread_mutex_lock(&write_mutex);
        assert(bgzf_write(struct_obj->bgzf_fp[0],struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
        if (PE==struct_obj->SeqType){
          assert(bgzf_write(struct_obj->bgzf_fp[1],struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);
        }
        pthread_mutex_unlock(&write_mutex);
        struct_obj->fqresult_r1->l =0;
        struct_obj->fqresult_r2->l =0;
      }
    }
    chr_idx = -1;
    iter++;
    localread++;
    current_reads_atom++;
  }
  if (struct_obj->bgzf_fp[0]){
    if (struct_obj->fqresult_r1->l > 0){
      pthread_mutex_lock(&write_mutex);
      assert(bgzf_write(struct_obj->bgzf_fp[0],struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
      if (PE==struct_obj->SeqType){assert(bgzf_write(struct_obj->bgzf_fp[1],struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
      pthread_mutex_unlock(&write_mutex);
      struct_obj->fqresult_r1->l =0;
      struct_obj->fqresult_r2->l =0;
    } 
  }

  for(int j=0; j<struct_obj->MaximumLength;j++){bam_destroy1(struct_obj->list_of_reads[j]);}
  
  free(sf->rand_alloc);
  delete sf;
  
  free(drand_alloc);
  free(drand_alloc_nt);
  free(drand_alloc_nt_adapt);
  free(drand_alloc_briggs);
  
  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);
  
  if(struct_obj->totalThreads>1)
    pthread_exit(NULL);
  
  return 0;
}
