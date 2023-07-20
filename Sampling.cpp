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
//#include <htslib/thread_pool.h>

#include <pthread.h>

#include "mrand.h"
#include "Briggs.h"
#include "Briggs2.h"
#include "NtSubModels.h"
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "ThreadGeneration.h"
#include "Sampling.h"
#include "sample_qscores.h"
#include "fasta_sampler.h"
#include "add_indels.h"
#include "NGSNGS_misc.h"

#define LENS 10000
#define MAXBINS 100
extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;


pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

void* Sampling_threads(void *arg) {
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  mrand_t *rand_alloc = mrand_alloc(struct_obj->rng_type,loc_seed);

  // sequence reads, original, modified, with adapters, pe
  char READ_ID[1024];
  char *FragmentSequence = (char*) calloc(LENS,1);
  char seq_r1[LENS] = {0};
  char seq_r2[LENS] = {0};
  char qual_r1[LENS] = "\0";
  char qual_r2[LENS] = "\0";
 
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
  
  //generating kstring for potential records of the stochastic indels
  char INDEL_INFO[1024];
  char INDEL_DUMP[1024];
  kstring_t *indel;
  indel =(kstring_t*) calloc(1,sizeof(kstring_t));
  indel->s = NULL;
  indel->l = indel->m = 0;

  size_t localread = 0;
  size_t current_reads_atom = 0;

  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc

  extern int SIG_COND;

  std::default_random_engine RndGen(loc_seed);
  
  sim_fragment *sf;
  if (struct_obj->LengthType==0)
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->FixedSize,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  else
    sf = sim_fragment_alloc(struct_obj->LengthType,struct_obj->distparam1,struct_obj->distparam2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  
  int C_to_T_counter = 0;int C_to_T_counter_rev = 0;
  int G_to_A_counter = 0;int G_to_A_counter_rev = 0;
  int refCp1 = 0;int refCTp1 = 0;int refCp2 = 0;int refCTp2 = 0;

  int modulovalue;
  if (reads > 1000000){
    modulovalue = 10;
  }
  else{ 
    modulovalue = 1;
  }
  size_t moduloread = reads/modulovalue;
  
  while (current_reads_atom < reads && SIG_COND){
    int posB = 0; int posE = 0;

    //sample fragmentlength
    int fraglength = getFragmentLength(sf); 

    if(fraglength < struct_obj->lowerlimit){
      fprintf(stderr,"lowerlimit%d\n",struct_obj->lowerlimit);
      fraglength = struct_obj->lowerlimit;
    }
    // Selecting genomic start position across the generated contiguous contigs for which to extract 
    int chr_idx = -1;

    //maxbases is the number of nucleotides we will work with.
    //set this to minimum of bases sequenced or fragment length.
    //if output is fasta then length is simply the entire fragment
    int maxbases = std::min(fraglength,struct_obj->maxreadlength);
    if (struct_obj->OutputFormat==faT ||struct_obj->OutputFormat==fagzT)
      maxbases = fraglength;

    //get shallow copy of chromosome, offset into, is defined by posB, and posE
    char *chrseq = sample(struct_obj->reffasta,rand_alloc,&chr,chr_idx,posB,posE,fraglength);
    //extracting a biological fragment of the reference genome
    memset(FragmentSequence,0,strlen(FragmentSequence));
    strncpy(FragmentSequence,chrseq+(posB),fraglength); // same orientation as reference genome 5' -------> FWD -------> 3'
    
    /*
    // not implemented yet - nor to be included within the near future as off december 2022
    if circ
      if (posE < posB+fraglenth) 
        memset(FragmentSequence,'0',strlen(FragmentSequence));
        strncpy(FragmentSequence,chrseq+(posB),fraglength);
        strncpy(FragmentSequence,chrseq,fraglength-(posE-posB));         
    */

    int skipread = 1;
    if(FragmentSequence[0]!='N' && FragmentSequence[(int)strlen(FragmentSequence)-1]!='N')
      skipread = 0;

    if(skipread==1){
      memset(FragmentSequence,0,strlen(FragmentSequence));
      continue;
    }

    int fragmentLength =posE-posB;
    assert(posE>=posB&&fragmentLength>20);
    //Generating random ID unique for each read output
    double rand_val_id = mrand_pop(rand_alloc);
    int rand_id = (rand_val_id * fraglength-1)+(rand_val_id*current_reads_atom);

    int strandR1 = mrand_pop(rand_alloc)>0.5?0:1;

    // Sequence alteration integers
    int FragMisMatch = 0;
    int has_seqerr = 0;    
    int has_indels = 0;

    // Stochastic structural variation model    
    if(struct_obj->DoIndel){
      double pars[4] = {struct_obj->IndelFuncParam[0],struct_obj->IndelFuncParam[1],struct_obj->IndelFuncParam[2],struct_obj->IndelFuncParam[3]};

      int ops[2] ={0,0};
      add_indel(rand_alloc,FragmentSequence,struct_obj->maxreadlength,pars,INDEL_INFO,ops);
      if (ops[0] > 0 && ops[1] == 0){
        has_indels = 1;
      }
      else if (ops[0] == 0 && ops[1] > 0){
        has_indels = 2;
      }
      else if (ops[0] > 0 && ops[1] > 0){
        has_indels = 3;
      }
      int IndelFragLen = strlen(FragmentSequence);
      maxbases = std::min(IndelFragLen,struct_obj->maxreadlength);
    }

    // Mismatch matrix input file
    if(struct_obj->doMisMatchErr){
      FragMisMatch = MisMatchFile(FragmentSequence,rand_alloc,struct_obj->MisMatch,struct_obj->MisLength);
    }

    if (strandR1 == 1){
      // 5' -------> REV -------> 3'
      ReversComplement(FragmentSequence);
    }

    char **FragRes;

    //Nucleotide alteration models only on the sequence itself which holds for fa,fq,sam
    int ReadDeam = 0;
    int Groupshift = mrand_pop(rand_alloc)>0.5?0:1; 
    int FragTotal = 4;
    int iter = 1; //iterating through all fragments
    if(struct_obj->DoBriggs){
      //For the none-biotin briggs model we need to store 4 fragments with slightly different deaminations patterns
      FragRes = new char *[FragTotal];
      for(int i=0;i<FragTotal;i++){
        FragRes[i] = new char[1024];
        memset(FragRes[i],'\0',1024);
      }
      
      ReadDeam=SimBriggsModel2(FragmentSequence, fragmentLength, 
        struct_obj->BriggsParam[0],
        struct_obj->BriggsParam[1],
        struct_obj->BriggsParam[2],
        struct_obj->BriggsParam[3],rand_alloc,FragRes,strandR1,
        C_to_T_counter,G_to_A_counter,C_to_T_counter_rev,G_to_A_counter_rev,
        refCp1,refCTp1,refCp2,refCTp2);

      if (struct_obj->Duplicates == 1){
        //keep one fragment out of the 4 possible
        Groupshift = mrand_pop_long(rand_alloc) % 4;
        FragTotal = Groupshift+1;
      }
      else if (struct_obj->Duplicates == 2){
        iter = 2;
      }
      else if (struct_obj->Duplicates == 4){
        // Keep all 4 fragments
        Groupshift = 0;
        FragTotal = 4;
      }
      
    }
    else if(struct_obj->DoBriggsBiotin){
      FragRes = new char *[1];
      FragRes[0] = FragmentSequence;
      Groupshift = 0;
      FragTotal = 1;
      ReadDeam=0;
      ReadDeam = SimBriggsModel(FragRes[0],fragmentLength,struct_obj->BriggsParam[0],
		    struct_obj->BriggsParam[1],
		    struct_obj->BriggsParam[2], 
		    struct_obj->BriggsParam[3],rand_alloc,
        strandR1,C_to_T_counter,G_to_A_counter,C_to_T_counter_rev,G_to_A_counter_rev);
    }
    else{
      // for the none-deaminated sequences we likewise need to generate PCR duplicates
      FragTotal = struct_obj->Duplicates;
      FragRes = new char *[FragTotal];
      for (int i = 0; i < FragTotal; i++){
        FragRes[i] = FragmentSequence;
      }
      if (struct_obj->Duplicates == 1){
        // keep one fragment out of the 4 possible
        FragTotal = 1;
      }
      else if (struct_obj->Duplicates == 2){
        //only if we want two duplicates we have to select a pair and iterate through the pair
        FragTotal = 2;
        Groupshift = 0;
      }
      else if (struct_obj->Duplicates == 4){
        //only if we want two duplicates we have to select a pair and iterate through the pair
        Groupshift = 0;
      } 
    }

    int chr_idx_array[struct_obj->Duplicates];
    for (int i = 0; i < struct_obj->Duplicates; i++){chr_idx_array[struct_obj->Duplicates]=chr_idx;}
    
    for (int FragNo = 0+Groupshift; FragNo < FragTotal; FragNo+=iter){
      qual_r1[0] = qual_r2[0] = seq_r1[0] = seq_r2[0] = '\0';

      //now copy the actual sequence into seq_r1 and seq_r2 if PE 
      strncpy(seq_r1,FragRes[FragNo],maxbases);

      if(PE==struct_obj->SeqType)
        strncpy(seq_r2,FragRes[FragNo]+(fraglength-maxbases),maxbases);
    
      //Remove reads which are all N, we remove both pairs if both reads are all N
  
      for(int i=0;skipread&&i<(int)strlen(seq_r1);i++)
        if(seq_r1[i]!='N')
          skipread = 0;
      
      if(PE==struct_obj->SeqType){
        for(int i=0; skipread && i<(int)strlen(seq_r2);i++)
          if(seq_r2[i]!='N')
          skipread = 0;
      }
      
      if(skipread==1){
        memset(qual_r1, 0, sizeof qual_r1); 
        memset(qual_r2, 0, sizeof qual_r2);
        memset(seq_r1, 0, sizeof seq_r1);
        memset(seq_r2, 0, sizeof seq_r2);
        continue;
      }
    
      if(strlen(seq_r1) < 20)
        continue;
      
      //generating sam output information
      int SamFlags[2] = {-1,-1}; //flag[0] is for read1, flag[1] is for read2
      if(struct_obj->DoBriggs){
        if (FragNo==0||FragNo==2){
          //The sequences are equal to the reference
          if (SE==struct_obj->SeqType){
            // R1  5' |---R1-->|--FWD------------> 3'
            //     3' |-----------REV------------> 3'
            SamFlags[0] = 0; // Forward strand
          }
          else if (PE==struct_obj->SeqType){
            // R1  5' |---R1-->|--FWD------------> 3'
            // R2  3' ------------REV---|<--R2---| 5'
            SamFlags[0] = 97; // Read paired, mate reverse strand, first in pair
            SamFlags[1] = 145; // Read paired, read reverse strand, second in pair
            ReversComplement(seq_r2);
          }
        }
        if (FragNo==1||FragNo==3){
          //The sequences are reverse complementary of the original reference orientation
          if (SE==struct_obj->SeqType){
            // R1  5' |-----------FWD------------> 3'
            //     3' |---R1-->|--REV------------> 3'
            SamFlags[0] = 16;
          }
          else if (PE==struct_obj->SeqType){
            //R2  5' ------------FWD---|<--R2---| 3'  
            //R1  3' |---R1-->|--REV------------> 5'
            SamFlags[0] = 81;
            SamFlags[1] = 161;
            ReversComplement(seq_r2);
          }
        }
      }
      else{
        if (SE==struct_obj->SeqType){
          if (strandR1 == 0)
            SamFlags[0] = 0;
          else if (strandR1 == 1){
            SamFlags[0] = 16;
          }
        }
        if (PE==struct_obj->SeqType){
          if (strandR1 == 0){
            SamFlags[0] = 97;
            SamFlags[1] = 145;
          }
          else if (strandR1 == 1){
            SamFlags[0] = 81;
            SamFlags[1] = 161;
          }
          ReversComplement(seq_r2);
        }
      }
      
      //so now everything is on 5->3 and some of them will be reverse complement to referene
      //now everything is the same strand as reference, which we call plus/+
    
      snprintf(READ_ID,1024,"T%d_RID%d_S%d_%s:%d-%d_length:%d_mod%d%d%d", struct_obj->threadno, rand_id,strandR1,chr,posB+1,posE,fraglength,ReadDeam,FragMisMatch,has_indels);

      int nsofts[2] = {0,0}; //this will contain the softclip information to be used by sam/bam/cram output format for adapter and polytail
      //below will contain the number of bases for R1 and R2 that should align to reference before adding adapters and polytail
      int naligned[2] = {(int)strlen(seq_r1),-1};
      if(PE==struct_obj->SeqType)
        naligned[1] = strlen(seq_r2);
      
      //add adapters
      if(struct_obj->AddAdapt){
        //Because i have reverse complemented the correct sequences and adapters depending on the strand origin (or flags), i know all adapters will be in 3' end
        nsofts[0] = std::min(struct_obj->maxreadlength-strlen(seq_r1),std::string(struct_obj->Adapter_1).length());
        strncpy(seq_r1+strlen(seq_r1),struct_obj->Adapter_1,nsofts[0]);
        /*
        Sampling.cpp: In function 'void* Sampling_threads(void*)':
        Sampling.cpp:369:16: warning: 'char* strncpy(char*, const char*, size_t)' specified bound depends on the length of the source argument [-Wstringop-truncation]
          369 |         strncpy(seq_r1+strlen(seq_r1),struct_obj->Adapter_1,nsofts[0]);
              |         ~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Sampling.cpp:368:77: note: length computed here
          368 |         nsofts[0] = std::min(struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
              |                                                                       ~~~~~~^~~~~~~~~~~~~~~~~~~~~~~
        Sampling.cpp:369:16: warning: 'char* strncpy(char*, const char*, size_t)' specified bound depends on the length of the source argument [-Wstringop-truncation]
          369 |         strncpy(seq_r1+strlen(seq_r1),struct_obj->Adapter_1,nsofts[0]);
              |         ~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Sampling.cpp:368:77: note: length computed here
          368 |         nsofts[0] = std::min(struct_obj->maxreadlength-strlen(seq_r1),strlen(struct_obj->Adapter_1));
        */
        if(PE==struct_obj->SeqType){
          nsofts[1] = std::min(struct_obj->maxreadlength-strlen(seq_r2),std::string(struct_obj->Adapter_2).length());
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
        if((int)strlen(seq_r1)!=(naligned[0]+nsofts[0])){
          fprintf(stderr,"Number of aligned bases + number of adap + poly does not match\n");
          exit(0);
        }
        //below only runs for PE that is when nalign[1] is not -1
        if(naligned[1]!=-1 && (int)strlen(seq_r2)!=naligned[1]+nsofts[1]){
          fprintf(stderr,"Number of aligned bases + number of adap + poly does not match\n");
          exit(0);
        }
      }

      //now seq_r1 and seq_r2 is completely done, we can generate the quality score if the format is different for Fasta
      if (struct_obj->OutputFormat==faT ||struct_obj->OutputFormat==fagzT){
        sprintf(READ_ID+strlen(READ_ID),"%d F%d",0,FragNo);
        ksprintf(fqs[0],">%s R1\n%s\n",READ_ID,seq_r1);//make this into read
        if (PE==struct_obj->SeqType)
        ksprintf(fqs[1],">%s R2\n%s\n",READ_ID,seq_r2);
      } 
      else{
        //Fastq and Sam needs quality scores
        has_seqerr = sample_qscores(seq_r1,qual_r1,strlen(seq_r1),struct_obj->QualDist_r1,struct_obj->NtQual_r1,rand_alloc,struct_obj->DoSeqErr,ErrProbTypeOffset);
        if (PE==struct_obj->SeqType)
          has_seqerr = sample_qscores(seq_r2,qual_r2,strlen(seq_r2),struct_obj->QualDist_r1,struct_obj->NtQual_r1,rand_alloc,struct_obj->DoSeqErr,ErrProbTypeOffset);
        
        sprintf(READ_ID+strlen(READ_ID),"%d F%d",has_seqerr,FragNo);

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

          if (struct_obj->Align){
            AlignCigar[0][0] = bam_cigar_gen(naligned[0], BAM_CMATCH);
            if(nsofts[0]>0){
              AlignCigar[0][1] = bam_cigar_gen(nsofts[0], BAM_CSOFT_CLIP);
              n_cigar[0] = 2;
            }
            if(SamFlags[0]==16||SamFlags[0]==81){
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
              if(SamFlags[0] == 97){
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
          insert = max_end - min_beg;

          const char* suffR1 = " R1";
          const char* suffR2 = " R2";

          char READ_ID2[1024];
          strcpy(READ_ID2,READ_ID);
          
          strcat(READ_ID,suffR1);
          int chr_max_end_mate = 0; //PNEXT 0-> unavailable for SE
          int insert_mate = 0; //TLEN
          int chr_idx_mate = -1;
          if(PE==struct_obj->SeqType){
            strcat(READ_ID2,suffR2);
            if (struct_obj->Align == 0){
              mapq = 255;
              SamFlags[0] = SamFlags[1] = 4;
              chr_idx_array[struct_obj->Duplicates] = -1;
              chr_max_end_mate = 0;
            }
            else{
              chr_max_end_mate = max_end;
              insert_mate = insert;
              chr_idx_mate = chr_idx_array[struct_obj->Duplicates];
            }        
          }
          if (struct_obj->Align == 0){
            mapq = 255;
            SamFlags[0] = SamFlags[1] = 4;
            chr_idx_array[struct_obj->Duplicates] = -1;
            min_beg = -1;
            max_end = 0;
          }

          //we have set the parameters accordingly above for no align and PE
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID),READ_ID,SamFlags[0],chr_idx_array[struct_obj->Duplicates],min_beg,mapq,
                n_cigar[0],AlignCigar[0],chr_idx_mate,chr_max_end_mate-1,insert_mate,strlen(seq_r1),seq_r1,qual_r1,l_aux);

          //write PE also
          if (PE==struct_obj->SeqType){
            bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],strlen(READ_ID2),READ_ID2,SamFlags[1],chr_idx_array[struct_obj->Duplicates],max_end-1,mapq,
              n_cigar[1],AlignCigar[1],chr_idx_array[struct_obj->Duplicates],min_beg,0-insert_mate,strlen(seq_r2),seq_r2,qual_r2,l_aux);
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
      
      if (struct_obj->DoIndel && struct_obj->IndelDumpFile != NULL){
        snprintf(INDEL_DUMP,1024,"%s\t%s\n",READ_ID,INDEL_INFO);
        ksprintf(indel,"%s",INDEL_DUMP);
        if (struct_obj->bgzf_fp[2]){
          if (indel->l > 0){
            pthread_mutex_lock(&write_mutex);
            assert(bgzf_write(struct_obj->bgzf_fp[2],indel->s,indel->l)!=0);
            pthread_mutex_unlock(&write_mutex);
            indel->l = 0;
          }
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
      memset(seq_r1, 0, sizeof seq_r1);
      memset(seq_r2, 0, sizeof seq_r2);

      localread++;
      current_reads_atom++;
      //printing out every tenth of the runtime
      if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
        fprintf(stderr,"\t-> Thread %d produced %zu reads with a current total of %zu\n",struct_obj->threadno,moduloread,current_reads_atom);
    }
    chr_idx = -1;
    
    if(struct_obj->DoBriggs){
      for(int i=0;i<4;i++)
        delete[] FragRes[i];
    }
    delete[] FragRes;
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
	    indel->l = 0;
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

  free(FragmentSequence);
  if(struct_obj->totalThreads>1)
    pthread_exit(NULL);
  
  return 0;
}