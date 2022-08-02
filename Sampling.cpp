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

#define LENS 4096
#define MAXBINS 100
unsigned char nuc2int[255];

pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

void* Sampling_threads(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Sampling_thread *struct_obj = (Parsarg_for_Sampling_thread*) arg;

  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno; 
  size_t genome_len = strlen(struct_obj->genome);
  mrand_t *drand_alloc = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_nt_adapt = mrand_alloc(struct_obj->rng_type,loc_seed);
  mrand_t *drand_alloc_briggs = mrand_alloc(struct_obj->rng_type,loc_seed);
  
  //if(fqT==struct_obj->OutputFormat|| fqgzT==struct_obj->OutputFormat)
  //  qualstringoffset = 33;
  
  // sequence reads, original, modified, with adapters, pe
  char seq_r1[1024] = {0};
  char seq_r1_mod[1024] = {0};

  char read[1024] = {0};
  char readadapt[1024] = {0};

  char seq_r2[1024] = {0};
  char seq_r2_mod[1024] = {0};
  char read2[1024] = {0};
  char readadapt2[1024] = {0};
  
  char Adapter_1[1024] = {0};
  char Adapter_2[1024] = {0};

  size_t reads = struct_obj -> reads;
  size_t BufferLength = struct_obj -> BufferLength;

  size_t rand_start;

  char qual_r1[1024] = "\0";
  char qual_r2[1024] = "\0"; // {0};

  std::string MonoPhosphateSeq;
  std::string MonoPhosphateQual;
  int MonoLen;
  if (struct_obj->PolyNt != 'F'){
    MonoPhosphateSeq = std::string(1000, struct_obj->PolyNt);
    MonoPhosphateQual = std::string(1000, '!'); //lowest nucleotide quality string
  }

  size_t localread = 0;
  int iter = 0;
  size_t current_reads_atom = 0;
  int readsizelimit;

  extern int SIG_COND;

  double distpar1 = struct_obj->distparam1;
  double distpar2 = struct_obj->distparam2;
  int LengthType = struct_obj->LengthType;

  std::default_random_engine RndGen(loc_seed);
  sim_fragment *sf;
  if (LengthType==0)
    sf = sim_fragment_alloc(LengthType,struct_obj->FixedSize,distpar2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  else
    sf = sim_fragment_alloc(LengthType,distpar1,distpar2,struct_obj->No_Len_Val,struct_obj->FragFreq,struct_obj->FragLen,struct_obj->rng_type,loc_seed,RndGen);
  
  size_t moduloread = reads/10;
  while (current_reads_atom < reads &&SIG_COND){
    if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
      fprintf(stderr,"\t-> Produced %zu reads with a current total of %zu\n",moduloread,current_reads_atom);
    // Fragment length generator 2000 
    int fraglength = getFragmentLength(sf);

    // Selecting genomic start position across the generated contiguous contigs for which to extract 
    double rand_val = mrand_pop(drand_alloc);
    rand_start = rand_val * (genome_len-300)+1; //genome_len-100000;
    int chr_idx = 0;
    if(struct_obj->chr_no>1){
      // identify which contig to sample from
      while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}//OBS MEMORY LEAK
      // Upper fragment length cap for those reads originating in close proximity to contig end position.
      if ((rand_start+fraglength)>struct_obj->size_cumm[chr_idx+1]){fraglength = struct_obj->size_cumm[chr_idx+1]-rand_start;} //OBS
    }

    //fprintf(stderr,"FRAGMENT LENGTH %d\n",fraglength);
    // Determing upper bound and its affect on readsize limit based on the fragment length¨

    //NB
    if(strcasecmp("false",struct_obj->QualFlag)==0){readsizelimit = fraglength;} //THIS WOULD BE THE CASE FOR FASATA
    else{readsizelimit = struct_obj->readcycle;}//OBS CHANGE SUCH THAT IT DOESNT HOLD TRUE FOR FASTA

    //Generating random ID unique for each read output
    double rand_val_id = mrand_pop(drand_alloc);
    int rand_id = (rand_val_id * fraglength-1); //100

    // extract the DNA sequence from the respective contig and start position when considering potential read size limitation given the read quality profile
    if (fraglength > readsizelimit)
      strncpy(seq_r1,struct_obj->genome+rand_start-1,readsizelimit);
    else
      strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);

    if(PE==struct_obj->SeqType){
      if (fraglength > readsizelimit)
	      strncpy(seq_r2,struct_obj->genome+rand_start+fraglength-1-readsizelimit,readsizelimit);
      else
	      strncpy(seq_r2,struct_obj->genome+rand_start-1,fraglength);
    }
    
    //Selecting strand, 0-> forward strand (+) 5'->3', 1 -> reverse strand (-) 3'->5'
    int strand = mrand_pop(drand_alloc)>0.5?0:1; //int strand = (int) (rand_start%2);
    
    //Remove reads which are all N? In this code we remove both pairs if both reads are all N
    int skipread = 1;
    for(int i=0;skipread&&i<(int)strlen(seq_r1);i++)
      if(seq_r1[i]!='N')
	      skipread = 0;
    for(int i=0;seq_r2 && skipread && i<(int)strlen(seq_r2);i++)
      if(seq_r2[i]!='N')
	      skipread = 0;
    
    // generating sam output information
    int seqlen = strlen(seq_r1);
    int flags[2] = {-1,-1}; //flag[0] is for read1, flag[1] is for read2
    
    // Immediately generating the proper orientation and corresponding flag for potential bam output
    if (SE==struct_obj->SeqType){
      if (strand == 0)
	      flags[0] = 0;
      else if (strand == 1){
        flags[0] = 16;
        ReversComplement(seq_r1);
      }
    }
    else if (PE==struct_obj->SeqType){
      if (strand == 0){
        flags[0] = 97;
        flags[1] = 145;
        ReversComplement(seq_r2);
        //fprintf(stderr,"STRAND %d AND FLAGS %d  %d \n",strand,flags[0],flags[1]);        
      }
      else if (strand == 1){
        //fprintf(stderr,"STRAND %d\n",strand);
        flags[0] = 81;
        flags[1] = 161;
        ReversComplement(seq_r1);
        //fprintf(stderr,"STRAND %d AND FLAGS %d  %d \n",strand,flags[0],flags[1]);        
      }
    }
    
    int Adapter1_len = 0;
    int Adapter2_len = 0;
    int SeqAdapt1_len = 0;
    int SeqAdapt2_len = 0;

    if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
      strncpy(Adapter_1, struct_obj->Adapter_1, sizeof(Adapter_1));//strncpy or memcpy
      Adapter1_len = strlen(Adapter_1);
      SeqAdapt1_len = seqlen+Adapter1_len;
      if (PE==struct_obj->SeqType){
        strncpy(Adapter_2, struct_obj->Adapter_2, sizeof(Adapter_2)); //strncpy or memcpy
        Adapter2_len = strlen(Adapter_2);
        //fprintf(stderr,"FIRST ADAPTER 2 LENGTH %d\n",Adapter2_len);
        SeqAdapt2_len = seqlen+Adapter2_len;
      }
      //for monophosphate the soft clip needs to be extended in length
      if (struct_obj->PolyNt != 'F'){
	      Adapter1_len = readsizelimit - seqlen; 
        Adapter2_len = readsizelimit - seqlen; 
        //fprintf(stderr,"FIRST ADAPTER 2 LENGTH WITH POLY %d\n",Adapter2_len);
      }
    }

    //Nucleotide alteration models only on the sequence itself which holds for fa,fq,sam
    if(strcasecmp(struct_obj->Briggs_flag,"true")==0){
      SimBriggsModel(seq_r1, seq_r1_mod, strlen(seq_r1),struct_obj->BriggsParam[0],
                                                        struct_obj->BriggsParam[1],
                                                        struct_obj->BriggsParam[2], 
                                                        struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
      strncpy(seq_r1, seq_r1_mod, sizeof(seq_r1));
        
      if (PE==struct_obj->SeqType){
        SimBriggsModel(seq_r2, seq_r2_mod, fraglength,struct_obj->BriggsParam[0], 
                                                      struct_obj->BriggsParam[1], 
                                                      struct_obj->BriggsParam[2], 
                                                      struct_obj->BriggsParam[3],loc_seed,drand_alloc_briggs);
        strncpy(seq_r2, seq_r2_mod, sizeof(seq_r2));
      }
    }
    
    if(strcasecmp("true",struct_obj->SubFlag)==0){
      MisMatchFile(seq_r1,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
      
      if (PE==struct_obj->SeqType)
	      MisMatchFile(seq_r2,drand_alloc_briggs,struct_obj->MisMatch,struct_obj->MisLength);
	    
    }
      
    char READ_ID[1024];int read_id_length;
    char *chrname = struct_obj->names[chr_idx];
      
      // Extract specific chromosome name similar to the bcf
    if (struct_obj->Variant_flag!=NULL && strcasecmp(struct_obj->Variant_flag ,"bcf")==0)
      chrname = struct_obj->names[0];
	    
    read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%zu-%zu_length:%d", struct_obj->threadno, rand_id,strand,chrname,
      rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],(int)fraglength);
    
    // Adding adapters before adding sequencing errors.
    if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
      // Because i have reverse complemented the correct sequences and adapters depending on the strand origin (or flags), i know all adapters will be in 3' end
      strcpy(read,seq_r1);
      strcat(read,Adapter_1);

      strcpy(read2,seq_r2);
      strcat(read2,Adapter_2);

      //saving both fasta and adapter to fasta format
      if (struct_obj->OutputFormat==faT ||struct_obj->OutputFormat==fagzT){
        ksprintf(struct_obj->fqresult_r1,">%s R1\n%s%s\n",READ_ID,seq_r1,Adapter_1);
        if (PE==struct_obj->SeqType)
	        ksprintf(struct_obj->fqresult_r2,">%s R2\n%s%s\n",READ_ID,seq_r2,Adapter_2);
	    }
      else{
        // Fastq and Sam needs quality score for both read and adapter, but the length cannot exceed the readcycle length inferred from read profile
        strncpy(readadapt, read, readsizelimit);
        if (PE==struct_obj->SeqType)
          strncpy(readadapt2, read2, readsizelimit);
          
        // Since the adapters don't have to be the same length, it is necessary to seperate Read 1 from Read 2 when generating quality string
        sample_qscores(readadapt,qual_r1,strlen(readadapt),struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt_adapt,struct_obj->ErrorFlag);

        if (struct_obj->PolyNt != 'F'){
          MonoLen = readsizelimit - strlen(readadapt);
	        if (struct_obj->OutputFormat==fqT ||struct_obj->OutputFormat==fqgzT)
            ksprintf(struct_obj->fqresult_r1,"@%s R1\n%s%s\n+\n%s%s\n",READ_ID,readadapt,MonoPhosphateSeq.substr(1,MonoLen).c_str(),qual_r1,MonoPhosphateQual.substr(1,MonoLen).c_str());
          else if (struct_obj->SAMout){
	          if(flags[0] == 0 || flags[0] == 97)
		        ksprintf(struct_obj->fqresult_r1,"%s%s",readadapt,MonoPhosphateSeq.substr(1,MonoLen).c_str());
            // CIGAR STRING YSXM, Y=ADAPTERLENGTH, X=SEQUENCE LENGTH, The sequence is equal to reference but adapter is reverse complementary
            else if(flags[0] == 16 || flags[0] == 81){     
              //change orientation back to reference genome
              ReversComplement(readadapt);
              //adding monophosphate
              ksprintf(struct_obj->fqresult_r1,"%.*s%.*s%s",seqlen,readadapt+Adapter1_len,Adapter1_len,readadapt,MonoPhosphateSeq.substr(1,MonoLen).c_str());
            }
          }
        }
        else{
          if (struct_obj->OutputFormat==fqT || struct_obj->OutputFormat==fqgzT)
            ksprintf(struct_obj->fqresult_r1,"@%s R1\n%s\n+\n%s\n",READ_ID,readadapt,qual_r1);
          else if (struct_obj->SAMout){
            if(flags[0] == 0 || flags[0] == 97)
              ksprintf(struct_obj->fqresult_r1,"%s",readadapt);
            else if(flags[0] == 16 || flags[0] == 81){
              ReversComplement(readadapt);
              ksprintf(struct_obj->fqresult_r1,"%.*s%.*s",seqlen,readadapt+Adapter1_len,Adapter1_len,readadapt);
            }
          }
        }

        // Adding sequencing errors for second read
        if (PE==struct_obj->SeqType){
          //generate nucleotide quality score
          sample_qscores(readadapt2,qual_r2,strlen(readadapt2),struct_obj->QualDist_r2,struct_obj->NtQual_r2,drand_alloc_nt_adapt,struct_obj->ErrorFlag);

          if (struct_obj->PolyNt != 'F'){
            MonoLen = readsizelimit - strlen(readadapt2);
            if (struct_obj->OutputFormat==fqT ||struct_obj->OutputFormat==fqgzT)
              ksprintf(struct_obj->fqresult_r2,"@%s R2\n%s%s\n+\n%s%s\n",READ_ID,readadapt2,MonoPhosphateSeq.substr(1,MonoLen).c_str(),qual_r2,MonoPhosphateQual.substr(1,MonoLen).c_str());
            else if (struct_obj->SAMout){
              if(flags[1] == 161)
                ksprintf(struct_obj->fqresult_r2,"%s%s",readadapt2,MonoPhosphateSeq.substr(1,MonoLen).c_str());
              else if(flags[1] == 145){
                ReversComplement(readadapt2);
                ksprintf(struct_obj->fqresult_r2,"%.*s%.*s%s",seqlen,readadapt2+Adapter2_len,Adapter2_len,readadapt2,MonoPhosphateSeq.substr(1,MonoLen).c_str());
              }
            }
          }
          else{
            if (struct_obj->OutputFormat==fqT ||struct_obj->OutputFormat==fqgzT)
              ksprintf(struct_obj->fqresult_r2,"@%s R2\n%s\n+\n%s\n",READ_ID,readadapt2,qual_r2);
            else if (struct_obj->SAMout){
              if(flags[1] == 161)
		            ksprintf(struct_obj->fqresult_r2,"%s",readadapt2);
              else if(flags[1] == 145){
                ReversComplement(readadapt2);
                ksprintf(struct_obj->fqresult_r2,"%.*s%.*s",seqlen,readadapt2+Adapter2_len,Adapter2_len,readadapt2);
              }
            }
          } 
        }
      }
    }
    // Saving reads without adapter
    else{
      //saving both fasta and adapter to fasta format
      if (struct_obj->OutputFormat==faT || struct_obj->OutputFormat==fagzT){
        ksprintf(struct_obj->fqresult_r1,">%s R1\n%s\n",READ_ID,seq_r1);
        if (PE==struct_obj->SeqType)
          ksprintf(struct_obj->fqresult_r2,">%s R2\n%s\n",READ_ID,seq_r2);
      }
      else{
        sample_qscores(seq_r1,qual_r1,seqlen,struct_obj->QualDist_r1,struct_obj->NtQual_r1,drand_alloc_nt,struct_obj->ErrorFlag);

        if (struct_obj->OutputFormat==fqT ||struct_obj->OutputFormat==fqgzT)
          ksprintf(struct_obj->fqresult_r1,"@%s R1\n%s\n+\n%s\n",READ_ID,seq_r1,qual_r1);
        else if (struct_obj->SAMout){
          if (flags[0] == 0 || flags[0] == 97)
	          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
          else if (flags[0] == 16 || flags[0] == 81){
            ReversComplement(seq_r1);
            ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
	        }
	      }

        if (PE==struct_obj->SeqType){
          sample_qscores(seq_r2,qual_r2,seqlen,struct_obj->QualDist_r2,struct_obj->NtQual_r2,drand_alloc_nt,struct_obj->ErrorFlag);

          if (struct_obj->OutputFormat==fqT ||struct_obj->OutputFormat==fqgzT)
            ksprintf(struct_obj->fqresult_r2,"@%s R2\n%s\n+\n%s\n",READ_ID,seq_r2,qual_r2);
          else if (struct_obj->SAMout){
            if (flags[1] == 161)
		          ksprintf(struct_obj->fqresult_r2,"%s",seq_r2);
            else if (flags[1] == 145){
		          ReversComplement(seq_r2);
              ksprintf(struct_obj->fqresult_r2,"%s",seq_r2);
	          }
          }
        }
      }
    }
    if (struct_obj->bgzf_fp[0]){
      if (struct_obj->fqresult_r1->l > BufferLength){
        pthread_mutex_lock(&write_mutex);
        assert(bgzf_write(struct_obj->bgzf_fp[0],struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
        if (PE==struct_obj->SeqType){assert(bgzf_write(struct_obj->bgzf_fp[1],struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
        pthread_mutex_unlock(&write_mutex);
        struct_obj->fqresult_r1->l =0;
        struct_obj->fqresult_r2->l =0;
      }
    }
    else if (struct_obj->SAMout){
      // Generate CI§GAR string for potential sam output with or without adapter
      //unaligned reads
      
      //Adapter1_len = strlen(Adapter_1);
      //SeqAdapt1_len = seqlen+Adapter1_len;
      uint32_t cigsunmap[2][10];
      size_t n_cigar_unmap=1;
      

      // Flags and cigar string for reads depending on adapters
      //prepare cigarstrings for 10 cigar operations
      uint32_t cigs[2][10];//cigs[0] is read1 cigs[1] is read2
      size_t n_cigar = 1;

      // By keeping the alignment information all reads will have matches in their cigar string     
      if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
        if(flags[0] == 0 || flags[0] == 97){
          // flag 0 is for SE so that would be same adapter (-a1) as provided for PE for flag 97
          if(readsizelimit>=SeqAdapt1_len){
            n_cigar = 2;
            cigs[0][0] = bam_cigar_gen(seqlen, BAM_CMATCH);
            cigs[0][1] = bam_cigar_gen(Adapter1_len, BAM_CSOFT_CLIP);
            
            cigsunmap[0][0] = bam_cigar_gen(SeqAdapt1_len, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
          }
          else if(readsizelimit>seqlen && readsizelimit<SeqAdapt1_len){
            // i.e fragment length of 148 and 10 adapter nucleotides will be large than inferred read
            n_cigar = 2;
            cigs[0][0] = bam_cigar_gen(seqlen, BAM_CMATCH);
            cigs[0][1] = bam_cigar_gen(Adapter1_len-(SeqAdapt1_len-readsizelimit), BAM_CSOFT_CLIP);
            //fprintf(stderr,"Read length %d and seq adapt length %d and soft clip %d\n",seqlen,SeqAdapt1_len,Adapter1_len-(SeqAdapt1_len-readsizelimit));
            cigsunmap[0][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
          }
          else{
            //reads above the inferred read length will simply be containing matches
            n_cigar = 1;
            cigs[0][0] = bam_cigar_gen(readsizelimit, BAM_CMATCH);
            cigsunmap[0][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
          }

          if(flags[1] == 145){ //strand == 1 for read 2
            if(readsizelimit>=SeqAdapt2_len){
              n_cigar = 2;
              cigs[1][1] = bam_cigar_gen(seqlen, BAM_CMATCH);
              cigs[1][0] = bam_cigar_gen(Adapter2_len, BAM_CSOFT_CLIP);
              cigsunmap[1][0] = bam_cigar_gen(SeqAdapt2_len, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
            }
            else if(readsizelimit>seqlen && readsizelimit<SeqAdapt2_len){
              // i.e fragment length of 148 and 10 adapter nucleotides will be large than inferred read
              n_cigar = 2;
              cigs[1][1] = bam_cigar_gen(seqlen, BAM_CMATCH);
              cigs[1][0] = bam_cigar_gen(Adapter2_len-(SeqAdapt2_len-readsizelimit), BAM_CSOFT_CLIP);
              //fprintf(stderr,"Read length %d and seq adapt length %d and soft clip %d\n",seqlen,SeqAdapt1_len,Adapter1_len-(SeqAdapt1_len-readsizelimit));
              cigsunmap[1][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
            }
            else{
              n_cigar = 1;
              cigs[1][0] = bam_cigar_gen(readsizelimit, BAM_CMATCH);
              cigsunmap[1][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
            }
          }
        }
        else if(flags[0] == 16 || flags[0] == 81){ //strand == 1 for read 1
          if(readsizelimit>=SeqAdapt1_len){
            n_cigar = 2;
            cigs[0][0] = bam_cigar_gen(Adapter1_len, BAM_CSOFT_CLIP);
            cigs[0][1] = bam_cigar_gen(seqlen, BAM_CMATCH);
            cigsunmap[0][0] = bam_cigar_gen(SeqAdapt1_len, BAM_CSOFT_CLIP);
          }
          else if(readsizelimit>seqlen && readsizelimit<SeqAdapt1_len){
            // i.e fragment length of 148 and 10 adapter nucleotides will be large than inferred read
            n_cigar = 2;
            cigs[0][0] = bam_cigar_gen(Adapter1_len-(SeqAdapt1_len-readsizelimit), BAM_CSOFT_CLIP);
            cigs[0][1] = bam_cigar_gen(seqlen, BAM_CMATCH);
            //fprintf(stderr,"Read length %d and seq adapt length %d and soft clip %d\n",seqlen,SeqAdapt1_len,Adapter1_len-(SeqAdapt1_len-readsizelimit));
            cigsunmap[0][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
          }
          else{
            n_cigar = 1;
            cigs[0][0] = bam_cigar_gen(readsizelimit, BAM_CMATCH);
            cigsunmap[0][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP);
          }

          if(flags[1] == 161){ //strand == 0 for read 2
            //fprintf(stderr,"Flag printf 3 %d & %d and strand %d\n",flags[0],flags[1],strand);
            if(readsizelimit>=SeqAdapt2_len){
              n_cigar = 2;
              cigs[1][0] = bam_cigar_gen(seqlen, BAM_CMATCH);
              cigs[1][1] = bam_cigar_gen(Adapter2_len, BAM_CSOFT_CLIP);
              cigsunmap[1][0] = bam_cigar_gen(SeqAdapt2_len, BAM_CSOFT_CLIP);
            }
            else if(readsizelimit>seqlen && readsizelimit<SeqAdapt2_len){
              // i.e fragment length of 148 and 10 adapter nucleotides will be large than inferred read
              n_cigar = 2;
              cigs[1][0] = bam_cigar_gen(seqlen, BAM_CMATCH);
              cigs[1][1] = bam_cigar_gen(Adapter2_len-(SeqAdapt2_len-readsizelimit), BAM_CSOFT_CLIP);
              //fprintf(stderr,"Read length %d and seq adapt length %d and soft clip %d\n",seqlen,SeqAdapt1_len,Adapter1_len-(SeqAdapt1_len-readsizelimit));
              cigsunmap[1][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP); //read 1 adapter 1 - no align
            }
            else{
              n_cigar = 1;
              cigs[1][0] = bam_cigar_gen(readsizelimit, BAM_CMATCH);
              cigsunmap[1][0] = bam_cigar_gen(readsizelimit, BAM_CSOFT_CLIP);
            }
          }
        }
      }
      else{
        n_cigar = 1;
        uint32_t cigar_bit_match_unmap  = bam_cigar_gen(seqlen, BAM_CSOFT_CLIP);
        uint32_t cigar_bit_match = bam_cigar_gen(seqlen, BAM_CMATCH);
        cigs[0][0] = cigar_bit_match;
        cigs[1][0] = cigar_bit_match;
        cigsunmap[0][0] = cigar_bit_match_unmap; // read 1 - no align
        cigsunmap[1][0] = cigar_bit_match_unmap; // read 2 - no align
        
      }
        //generating id, position and the remaining sam field information
      size_t l_aux = 2; uint8_t mapq = 60;
      hts_pos_t min_beg, max_end, insert; //max_end, insert;
      min_beg = rand_start-struct_obj->size_cumm[chr_idx] - 1;
      max_end = rand_start-struct_obj->size_cumm[chr_idx] + fraglength - 1;
      insert = max_end - min_beg + 1;

      // Change the orientation of quality string similar to changing read orientation
      if (flags[0] == 16 || flags[0] == 81)
	      reverseChar(qual_r1,strlen(seq_r1));
      else if (flags[1] == 145)
	      reverseChar(qual_r2,strlen(seq_r2));
        
      char READIDR1[1024];
	    char READIDR2[1024];
      const char* suffR1 = " R1";const char* suffR2 = " R2";
      strcpy(READIDR1,READ_ID);strcat(READIDR1,suffR1);

      if (SE==struct_obj->SeqType){
        //Utlizing the sam format as a sequence container
        if (struct_obj->NoAlign == 'T'){
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR1),READIDR1,4,-1,-1,
          255,n_cigar_unmap,cigsunmap[0],-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,NULL,0);
        }
        else{ //saving 'alignment' information
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR1),READIDR1,flags[0],chr_idx,min_beg,mapq,
          n_cigar,cigs[0],-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
        }
        //const char* MDtag = "\tMD:Z:"; //bam_aux_update_str(struct_obj->list_of_reads[struct_obj->LengthData-1],"MD",5,"MDtag");
      }
      else if (PE==struct_obj->SeqType){
        strcpy(READIDR2,READ_ID);strcat(READIDR2,suffR2);
        if (struct_obj->NoAlign == 'T'){
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR1),READIDR1,4,-1,-1,
          255,n_cigar_unmap,cigsunmap[0],-1,-1,0,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,NULL,0);
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR2),READIDR2,4,-1,-1,
          255,n_cigar_unmap,cigsunmap[1],-1,-1,0,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,NULL,0);
        }
        else{
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR1),READIDR1,flags[0],chr_idx,min_beg,mapq,
          n_cigar,cigs[0],chr_idx,max_end,insert,strlen(struct_obj->fqresult_r1->s),struct_obj->fqresult_r1->s,qual_r1,l_aux);
	    
          bam_set1(struct_obj->list_of_reads[struct_obj->LengthData++],read_id_length+strlen(suffR2),READIDR2,flags[1],chr_idx,max_end,mapq,
          n_cigar,cigs[1],chr_idx,min_beg,0-insert,strlen(struct_obj->fqresult_r2->s),struct_obj->fqresult_r2->s,qual_r2,l_aux);
        }
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
    memset(qual_r1, 0, sizeof qual_r1); 
    memset(qual_r2, 0, sizeof qual_r2);  
    memset(seq_r1, 0, sizeof seq_r1);
    memset(seq_r1_mod, 0, sizeof seq_r1_mod);

    memset(seq_r2, 0, sizeof seq_r2);
    memset(seq_r2_mod, 0, sizeof seq_r2_mod);
      
    memset(readadapt, 0, sizeof readadapt);
    memset(readadapt2, 0, sizeof readadapt2);

      // didnt work memset(Adapter_1, 0, sizeof Adapter_1);memset(Adapter_2, 0, sizeof Adapter_2);
            
    chr_idx = 0;
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

  free(struct_obj->size_cumm);
  free(struct_obj->names);
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