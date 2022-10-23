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
#include "add_indels.h"
#include "NGSNGS_misc.h"

#define LENS 4096
#define MAXBINS 100
extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;


pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  mrand_t *rand_alloc = mrand_alloc(struct_obj->rng_type,loc_seed);

  char *seq;//actual sequence, this is unallocated
  int fraglength;

  // sequence reads, original, modified, with adapters, pe
  char READ_ID[1024];
  char FragmentSequence[4096];
  char seq_r1[1024] = {0};
  char seq_r2[1024] = {0};
  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};
 
  
  size_t reads = struct_obj -> reads;
  //fprintf(stderr,"INSIDE EACH THREAD NUMBER OF READS %zu \n",reads);
  size_t BufferLength = struct_obj -> BufferLength;

  int ErrProbTypeOffset = 0;
  if (struct_obj->OutputFormat==fqT || struct_obj->OutputFormat==fqgzT){ErrProbTypeOffset=33;}
  kstring_t *fqs[2];
  for(int i=0;i<2;i++){
    fqs[i] =(kstring_t*) calloc(1,sizeof(kstring_t));
    fqs[i]->s = NULL;
    fqs[i]->l = fqs[i]->m = 0;
  }
  
  //generating kstring for potential records of the stochastic indels
  char INDEL_INFO[512];
  char INDEL_DUMP[1024];
  kstring_t *indel;
  indel =(kstring_t*) calloc(1,sizeof(kstring_t));
  indel->s = NULL;
  indel->l = indel->m = 0;

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
  
  int C_total = 0;int C_to_T_counter = 0;int C_to_T_counter_rev = 0;int C_total_rev=0;
  int G_total = 0;int G_to_A_counter = 0;int G_to_A_counter_rev = 0;int G_total_rev=0;
  
  int modulovalue;
  if (reads > 1000000){
    modulovalue = 10;
  }
  else{ 
    modulovalue = 1;
  }
  size_t moduloread = reads/modulovalue;
  
  while (current_reads_atom < reads && SIG_COND){
    //lets start by resetting out datastructures to NULL, nill, nothing.
    qual_r1[0] = qual_r2[0] = seq_r1[0] = seq_r2[0] = '\0'; //Disse skal jo rykkes hvis vi bruger et char** til fragmenter
    int posB = 0; int posE = 0;//this is the first and last position of our fragment

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
    double rand_val_id = mrand_pop(rand_alloc);//<- SHIT
    //fprintf(stderr,"RAndom values for ID %lf \t %lf\n",rand_val_id,mrand_pop(drand_alloc));
    int rand_id = (rand_val_id * fraglength-1)+(rand_val_id*current_reads_atom); //100
    //why above? just to get random unique?

    int strandR1 = mrand_pop(rand_alloc)>0.5?0:1; //int strand = (int) (rand_start%2);

    //get shallow copy of chromosome, offset into, is defined by posB, and posE
    //fprintf(stderr,"NEW READ \n");
    //fprintf(stderr,"beg %d end %d fragleng %d max %d \n",posB,posE,fraglength,maxbases);
    char *chrseq = sample(struct_obj->reffasta,rand_alloc,&chr,chr_idx,posB,posE,fraglength);
    //extracting a biological fragment of the reference genome
    //fprintf(stderr,"beg %d end %d len: %lu frag %lu \n",posB,posE,strlen(chrseq),fraglength);

    strncpy(FragmentSequence,chrseq+(posB),(posE-posB)); // same orientation as reference genome 5' -------> FWD -------> 3'
    int fragmentLength = strlen(FragmentSequence);
    
    // Sequence alteration integers
    int FragDeam = 0;
    int has_indels = 0;
    int FragMisMatch = 0;

    // Nucleotide alteration models

    // Stochastic structural variation model    
    if(struct_obj->DoIndel){
      double pars[4] = {struct_obj->IndelFuncParam[0],struct_obj->IndelFuncParam[1],struct_obj->IndelFuncParam[2],struct_obj->IndelFuncParam[3]};
      //fprintf(stderr,"DO INDELS");
      /*if (struct_obj->DoIndel && struct_obj->IndelDumpFile != NULL){
        fprintf(stderr,"FILE DUMP TEST %s \n",struct_obj->IndelDumpFile);
      }*/
      //fprintf(stderr,"BEFORE FUNCTION %s\n",INDEL_INFO);
      fprintf(stderr,"\n1) %s len: %lu\n",FragmentSequence,strlen(FragmentSequence));
      has_indels = add_indel(rand_alloc,FragmentSequence,struct_obj->maxreadlength,pars,INDEL_INFO);
      fprintf(stderr,"\n2) %s len: %lu\n",FragmentSequence,strlen(FragmentSequence));
      int IndelFragLen = strlen(FragmentSequence);
      maxbases = std::min(IndelFragLen,struct_obj->maxreadlength);
      //fprintf(stderr,"AFTER FUNCTION\n%s",INDEL_INFO);
    }

    // Mismatch matrix input file
    if(struct_obj->doMisMatchErr){
      FragMisMatch = MisMatchFile(FragmentSequence,rand_alloc,struct_obj->MisMatch,struct_obj->MisLength);
    }
    //./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 10 -s 1 -l 2000 -seq SE -indel 0.05,0.1,0.1,0.2 -q1 Test_Examples/AccFreqL150R1.txt -f fq -o MycoBactBamSEOut
    //./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 1 -s 1 -l 100 -seq SE -indel 0.05,0.1,0.1,0.2 -q1 Test_Examples/AccFreqL150R1.txt -f fq -o MycoBactBamSEOut
    //generates insertions and deletions to the original fragment, before extracting R1 and R2
    
    if (strandR1 == 1){
      // 5' -------> REV -------> 3'
      ReversComplement(FragmentSequence);
    }

    /*

    INSERT NYE BRIGGS HER, LAV CHAR** OGSÅ FOR GAMLE BRIGGS OG SÅ GEM LÆNGDEN OG SÅ LAV ET TIL
    SVARENDE FOR LOOP SOM I DIN GAMLE SAMPLING OG LØB IGENNEM DEM ALLESAMMEN
    */


    //fprintf(stderr,"chr %s \t idx %d \n",chr,chr_idx);
    //now copy the actual sequence into seq_r1 and seq_r2 if PE 
    strncpy(seq_r1,FragmentSequence,maxbases);
    //fprintf(stderr,"sequence %s\n",seq_r1);
    //simulate indels for the fragment

    if(PE==struct_obj->SeqType)
      strncpy(seq_r2,FragmentSequence+(fraglength-maxbases),maxbases);
      //fprintf(stderr,"sequence pos is %d \t %s\n",strlen(FragmentSequence)-maxbases,seq_r2);

    //fprintf(stderr,"CHECKING THE LENGTH ISSUES %d \t %d \t %d \t %d \n",fraglength,struct_obj->maxreadlength,maxbases,strlen(seq_r1));

    /*
      ||------R1------>,,,,,,,,,,|-------R2----->||
    */
    
    //fprintf(stderr,"Current read %zu and max bases %d and fraglent %d  and coordinates %d \t %d \t seq length %d\n",current_reads_atom,maxbases,fraglength,posB,posE,strlen(seq_r1));
    //Selecting strand, 0-> forward strand (+) 5'->3', 1 -> reverse strand (-) 3'->5'
    //rename to strand to strandR1
    
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
    // NB SKAL RYKKES NED TIL AT PASSE MED DE NYE CHAR ARRAYS EFTER BRIGGS ER LAVET
    
    // generating sam output information
    int seqlen = strlen(seq_r1);

    if(strlen(seq_r1) < 20)
      continue;

    int SamFlags[2] = {-1,-1}; //flag[0] is for read1, flag[1] is for read2
    
    //now everything is the same strand as reference, which we call plus/+
    //lets flip to 5 to 3
    if (SE==struct_obj->SeqType){
      if (strandR1 == 0){
	      SamFlags[0] = 0;
        if (seq_r1[0]=='C'||seq_r1[0]=='c'){C_total++;}
        if (seq_r1[seqlen-1]=='G'||seq_r1[seqlen-1]=='g'){G_total++;}
      }
      else if (strandR1 == 1){
        SamFlags[0] = 16;
        ReversComplement(seq_r1);
        if (seq_r1[0]=='C'||seq_r1[0]=='c'){C_total_rev++;}
        if (seq_r1[seqlen-1]=='G'||seq_r1[seqlen-1]=='g'){G_total_rev++;}
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
    int ReadDeam = 0;
    if(struct_obj->DoBriggs || struct_obj->DoBriggsBiotin){
      /*fprintf(stderr,"-----------\nDO Briggs model %d \n",struct_obj->DoBriggs);
      fprintf(stderr,"-----------\nDO Briggs Biotin model %d \n",struct_obj->DoBriggsBiotin);
      fprintf(stderr,"BRIGGS PARAMETERS %f\n",struct_obj->BriggsParam[0]);*/
      ReadDeam=0;
      ReadDeam = SimBriggsModel(seq_r1,fraglength,struct_obj->BriggsParam[0],
		     struct_obj->BriggsParam[1],
		     struct_obj->BriggsParam[2], 
		     struct_obj->BriggsParam[3],rand_alloc,strandR1,C_to_T_counter,G_to_A_counter,C_to_T_counter_rev,G_to_A_counter_rev);
      
      if (struct_obj->DoBriggs){
        int RevBriggs = mrand_pop(rand_alloc)>0.5?0:1;
        //fprintf(stderr,"INSIDE DO BRIGGS LOOP %d\n",RevBriggs);
        // we consider the different strands using the meyer 2010 article in 50% of the cases
        if(RevBriggs == 1){
          //fprintf(stderr,"INSIDE REVBRIGGS IF with strand %d\n",strandR1);
          if (strandR1 == 0){
            if(struct_obj->SAMout){
              //fprintf(stderr,"INSIDE SAM FOR STRAND 0 IF\n");
              //in the sam output we shouldn't change the direction to be 3'->5', so we only change the flag
              SamFlags[0] = 16;
            }
            else{
              //fprintf(stderr,"INSIDE FQ FOR STRAND 0 IF\n");
              // Here we need to change the orientation
              ReversComplement(seq_r1);
            }
          }
          else if (strandR1 == 1){
            if(struct_obj->SAMout){
              SamFlags[0] = 0;
            }
            else{
              ReversComplement(seq_r1);
            }         
          }
        }
        //otherwise we continue with the same strand orientation such that the forward reads aren't affected, and the reverse reads are only reversecomplemented for the sam output structure
      }
      
      if (PE==struct_obj->SeqType){
        ReadDeam = SimBriggsModel(seq_r2, fraglength,struct_obj->BriggsParam[0], 
		       struct_obj->BriggsParam[1], 
		       struct_obj->BriggsParam[2], 
		       struct_obj->BriggsParam[3],rand_alloc,strandR1,C_to_T_counter,G_to_A_counter,C_to_T_counter_rev,G_to_A_counter_rev);
      }
    } 
      
    snprintf(READ_ID,1024,"T%d_RID%d_S%d_%s:%d-%d_length:%d_mod%d%d%d%d", struct_obj->threadno, rand_id,strandR1,chr,posB+1,posE,fraglength,ReadDeam,FragMisMatch,0,0);

    if (struct_obj->DoIndel && struct_obj->IndelDumpFile != NULL){
      //ksprintf(indel,"%s\t%s\n",READ_ID,INDEL_INFO);
      snprintf(INDEL_DUMP,1024,"%s\t%s\n",READ_ID,INDEL_INFO);
      //fprintf(stderr,"INFO INDEL ID %s",INDEL_DUMP);
      ksprintf(indel,"%s",INDEL_DUMP);
      if (struct_obj->bgzf_fp[2]){
        if (indel->l > 0){
          pthread_mutex_lock(&write_mutex);
          assert(bgzf_write(struct_obj->bgzf_fp[2],indel->s,indel->l)!=0);
          pthread_mutex_unlock(&write_mutex);
          indel->l = indel->l = 0;
        }
      }
    }
    
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
      sample_qscores(seq_r1,qual_r1,strlen(seq_r1),struct_obj->QualDist_r1,struct_obj->NtQual_r1,rand_alloc,struct_obj->DoSeqErr,ErrProbTypeOffset);
      //if(strandR1==0){fprintf(stderr,"%s\n",seq_r1);}
      if (PE==struct_obj->SeqType)
      	sample_qscores(seq_r2,qual_r2,strlen(seq_r2),struct_obj->QualDist_r1,struct_obj->NtQual_r1,rand_alloc,struct_obj->DoSeqErr,ErrProbTypeOffset);
      
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
      insert = max_end - min_beg;
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
      int chr_idx_mate = -1; // RNEXT 0 -> "=" -1 -> "*s"
      int chr_max_end_mate = 0; //PNEXT 0-> unavailable for SE
      int insert_mate = 0; //TLEN
      if(PE==struct_obj->SeqType){
        strcat(READ_ID2,suffR2);
        if (struct_obj->NoAlign == 0){
          mapq = 255;
          SamFlags[0] = SamFlags[1] = 4;
          chr_idx = -1;
          chr_max_end_mate = 0;
        }
        else{
          chr_idx_mate = chr_idx;
          chr_max_end_mate = max_end;
          insert_mate = insert;
          //fprintf(stderr,"CHR IDX %d\n",chr_idx_mate);
        }        
      }
      if (struct_obj->NoAlign == 0){
        mapq = 255;
        SamFlags[0] = SamFlags[1] = 4;
        chr_idx = -1;
        min_beg = -1;
        max_end = 0;
      }

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
        }
        pthread_mutex_unlock(&write_mutex);
        struct_obj->LengthData = 0;
      }
      fqs[0]->l = fqs[1]->l = 0;
      }
      
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
    memset(FragmentSequence,0,sizeof FragmentSequence);  
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

  if (struct_obj->bgzf_fp[2]){
    if (indel->l > 0){
      pthread_mutex_lock(&write_mutex);
      assert(bgzf_write(struct_obj->bgzf_fp[0],indel->s,indel->l)!=0);
      pthread_mutex_unlock(&write_mutex);
	    indel->l = indel->l = 0;
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
  
  free(indel->s);
  free(indel);

  free(rand_alloc);
  
  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);

  if(struct_obj->DoBriggs || struct_obj->DoBriggsBiotin){
    fprintf(stdout,"C>T and G>A frequency at position 1 5' and 1 3' for forward strand\t%f\t%f\n",(double)C_to_T_counter/(double)C_total,(double)G_to_A_counter/(double)G_total);
    fprintf(stdout,"C>T and G>A frequency at position 1 5' and 1 3' for reverse strand\t%f\t%f\n",(double)C_to_T_counter_rev/(double)C_total_rev,(double)G_to_A_counter_rev/(double)G_total_rev);
  }

  if(struct_obj->totalThreads>1)
    pthread_exit(NULL);
  
  return 0;
}
