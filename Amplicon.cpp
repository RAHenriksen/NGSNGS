#include "Briggs.h"
#include "mrand.h"
#include "fasta_sampler.h"
#include "NtSubModels.h"
#include "add_indels.h"
#include "HelpPage.h"
#include "Amplicon_cli.h"
#include "Amplicon.h"
#include "NGSNGS_misc.h"
#include "version.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <iostream>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

#define LENS 4096

pthread_mutex_t amplicon_write_mutex = PTHREAD_MUTEX_INITIALIZER;

KSEQ_INIT(BGZF*, bgzf_read);


void* AmpliconThreadInitialize(ampliconformat_e OutFormat,const char* Amplicon_in_fp,
  int filetype,const char* Amplicon_out_fp,
  const char* Subprofile,int threads,
  char *Briggs,char *Indel,int seed,
  int rng_type,size_t moduloread,size_t totalLines,size_t linesPerThread) {
  
  // initialize file type for output name to convert files
  const char* suffix = NULL;
  const char* mode = NULL;

  char fileout[512];
  const char* fileprefix = Amplicon_out_fp;
  strcpy(fileout,fileprefix);

  switch(OutFormat){
    case faT:
      mode = "wu";
      suffix = ".fa";
      break;
    case fagzT:
      mode = "wb";
      suffix = ".fa.gz";
      break;
    case fqT:
      mode = "wu";
      suffix = ".fq";
      break;
    case fqgzT:
      mode = "w";
      suffix = ".fq.gz";
      break;
    case samT:
      mode = "ws";
      suffix = ".sam";
      break;
    case bamT:
      mode = "wb";
      suffix = ".bam";
      break;
    default:
      fprintf(stderr,"\t-> Fileformat is currently not supported \n");
      break;
  }

  strcat(fileout,suffix);

  //fprintf(stderr,"\t-> File output name is %s\n",fileout);
  const char* filename_out = fileout;
  
  //generating mismatch matrix to parse for each string
  double* MisMatchFreqArray = new double[LENS];
  int mismatchcyclelength = 0;
  int doMisMatchErr = 0;
  if (Subprofile != NULL){
    MisMatchFreqArray = MisMatchFileArray(MisMatchFreqArray,Subprofile,mismatchcyclelength);
    doMisMatchErr = 1;
  }


  // overvej her at definer dine BGZF som null og så afhængigt af parametrene så kan du åbne den, også definerer den? altså så kan du netop også nemt gøre det med bam, ellers så bare lave det nemt og flot med et seperate script
  BGZF* amplicon_in = NULL;
  BGZF* amplicon_out = NULL;
  samFile *amplicon_in_sam = NULL;
  samFile *amplicon_out_sam = NULL;
  bam_hdr_t *hdr = NULL;
  
  if(filetype<2){
    // input fa or fastq
    amplicon_in = bgzf_open(Amplicon_in_fp, "r");  
  }
  else{
    //input sequence alignment/map format
    amplicon_in_sam = hts_open(Amplicon_in_fp, "r");

    // Read the header
    hdr = sam_hdr_read(amplicon_in_sam);
    //bam_header_t *sam_header_read(
    if (hdr == NULL) {
      fprintf(stderr, "Error reading header from %s\n", Amplicon_in_fp);
      sam_close(amplicon_in_sam);
      exit(1);
    }
  }

  if(amplicon_in == NULL && amplicon_in_sam == NULL){
    fprintf(stderr, "Error reading header from %s\n", Amplicon_in_fp);
    exit(1);
  }

  if(OutFormat == samT || OutFormat == bamT){
    //output file
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
    //const char* output_bam = "test.bam";
    amplicon_out_sam = sam_open_format(filename_out, mode, fmt_hts);
    if (sam_hdr_write(amplicon_out_sam, hdr) < 0) {
      fprintf(stderr, "Error writing header to BAM file %s\n", filename_out);
      bam_hdr_destroy(hdr);
      sam_close(amplicon_out_sam);
      exit(1);
    }
  }
  else{
    amplicon_out = bgzf_open(filename_out, mode);
    if (amplicon_out == NULL) {
      fprintf(stderr, "Error opening output file.\n");
      exit(1);
    }
  }

  // Create an array to hold thread IDs dynamically
  pthread_t* mythreads = new pthread_t[threads];

  // Create an array to hold thread arguments dynamically
  struct_for_amplicon_threads* amp_thread_struct = new struct_for_amplicon_threads[threads];

  // Create threads dynamically
  for (int i = 0; i < threads; ++i) {
    int startLine = i * linesPerThread+1;
    int endLine = (i == threads - 1) ? totalLines : startLine + linesPerThread - 1;

    //fprintf(stderr, "Thread %d reads in line start %d to line end %d\n", i, startLine, endLine);
    amp_thread_struct[i].OutFormat = OutFormat;
    amp_thread_struct[i].filetype = filetype;
    amp_thread_struct[i].moduloread = moduloread;
    amp_thread_struct[i].totalLines = totalLines;
    amp_thread_struct[i].linesPerThread = linesPerThread;

    amp_thread_struct[i].amplicon_in_fp = amplicon_in;
    amp_thread_struct[i].amplicon_out_fp = amplicon_out;

    amp_thread_struct[i].amplicon_in_sam = amplicon_in_sam;
    amp_thread_struct[i].amplicon_out_sam = amplicon_out_sam;
    amp_thread_struct[i].hdr = hdr;

    amp_thread_struct[i].Briggs = Briggs;
    amp_thread_struct[i].startLine = startLine;
    amp_thread_struct[i].endLine = endLine;
    amp_thread_struct[i].threadid = i;
    amp_thread_struct[i].Seed = seed; 
    amp_thread_struct[i].rng_type = rng_type;

    // 3) misincorporation matrix
    amp_thread_struct[i].MisMatch = MisMatchFreqArray;
    amp_thread_struct[i].MisLength = (int) mismatchcyclelength;
    amp_thread_struct[i].doMisMatchErr = doMisMatchErr;

    // Random indels
    amp_thread_struct[i].Indel = Indel;

    // Create each thread
     
    if(filetype<2){
      pthread_create(&mythreads[i], NULL, ProcessFAFQ, &amp_thread_struct[i]);
    }
    else{
      pthread_create(&mythreads[i], NULL, ProcessBAM, &amp_thread_struct[i]);
    }
    
  }
  
  // Wait for all threads to finish
  for (int i = 0; i < threads; ++i) {
    pthread_join(mythreads[i], NULL);
    //fprintf(stderr, "Thread %d finished\n", i+1);
  }

  // Delete allocated memory
  delete[] mythreads;
  delete[] amp_thread_struct;


  if(OutFormat == samT || OutFormat == bamT){
    sam_close(amplicon_out_sam);
  }
  else{
    // Close the files
    bgzf_close(amplicon_out);
  }

  if(filetype<2){
    // input fa or fastq
    bgzf_close(amplicon_in);  
  }
  else{
    sam_hdr_destroy(hdr);
    sam_close(amplicon_in_sam);
  }
  return NULL;
}


void* ProcessFAFQ(void* args){
  struct_for_amplicon_threads* amp_thread_struct = (struct_for_amplicon_threads*)args;

  BGZF* amplicon_in_fp = amp_thread_struct->amplicon_in_fp;
  BGZF* amplicon_out_fp = amp_thread_struct->amplicon_out_fp;
  int startLine = amp_thread_struct->startLine;
  int endLine = amp_thread_struct->endLine;
  int threadid = amp_thread_struct->threadid;

  int filetype = amp_thread_struct->filetype;
  
  //char* Briggs = amp_thread_struct->Briggs;
  //char* IndelInputParam = amp_thread_struct->Indel;

  long int Seed = amp_thread_struct->Seed;

  //fprintf(stderr,"initialize thread %d reading chunk starting from line %d to ending line %d\n",threadid,startLine,endLine);
  
  // Initialize briggs parameters
  float Param[4];
  if (amp_thread_struct->Briggs != NULL){
    char* BriggsParam;
    BriggsParam = strdup(amp_thread_struct->Briggs);
    Param[0] = atof(strtok(BriggsParam,"\", \t"));
    Param[1] = atof(strtok(NULL,"\", \t"));
    Param[2] = atof(strtok(NULL,"\", \t"));
    Param[3] = atof(strtok(NULL,"\", \t"));
      
    free(BriggsParam);
  }

  // initialize indel parameters
  float IndelFuncParam[4];
  if (amp_thread_struct->Indel != NULL){
    char* IndelInputParam = strdup(amp_thread_struct->Indel);
    IndelFuncParam[0] = atof(strtok(IndelInputParam,"\", \t"));
    IndelFuncParam[1] = atof(strtok(NULL,"\", \t"));
    IndelFuncParam[2] = atof(strtok(NULL,"\", \t"));
    IndelFuncParam[3] = atof(strtok(NULL,"\", \t"));
    
    free(IndelInputParam); 
  }

  //generating kstring for potential records of the stochastic indels
  char INDEL_INFO[1024];  

  // Count the number of processed reads  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);
  size_t moduloread = amp_thread_struct->moduloread;
  size_t localread = 0;
  size_t current_reads_atom = 0;

  //allocate the random generator
  mrand_t *mr = mrand_alloc(amp_thread_struct->rng_type,Seed);

  // Seek to the appropriate starting line
  kseq_t *FQseq = kseq_init(amplicon_in_fp);

  int ErrProbTypeOffset=33; //quality score offset, 33 only for fastq it should be 0 for bam

  while (startLine <= endLine){
    //fprintf(stderr,"line %d\n",startLine);
    if (kseq_read(FQseq) < 0) {
      //fprintf(stderr, "Error reading sequence in thread %d\n", threadid);
      break;    
    }

    localread++;
    current_reads_atom++;
    //printing out every tenth of the runtime
    if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
      fprintf(stderr,"\t-> Processed %zu reads with a current total of %zu\n",moduloread,current_reads_atom);

    // Sequence alteration integers
    int FragMisMatch = 0;
    int has_indels = 0;
    int ReadDeam=0;

    // deamination
    if (amp_thread_struct->Briggs != NULL){
      int strand = mrand_pop(mr)>0.5?0:1;
      ReadDeam = SimBriggsModel_amplicon(&FQseq->seq,Param[0],Param[1],Param[2],Param[3], mr);
    }

    // Mismatch matrix input file
    if(amp_thread_struct->doMisMatchErr > 0){
      //fprintf(stderr,"INSIDE mf\n");    
      FragMisMatch = MisMatchFile_kstring(&FQseq->seq,mr,amp_thread_struct->MisMatch,amp_thread_struct->MisLength);
      //fprintf(stderr,"FragMisMatch val %d \n",FragMisMatch);
    }

    // Stochastic structural variation model    
    if(amp_thread_struct->Indel != NULL){
      double pars[4] = {IndelFuncParam[0],IndelFuncParam[1],IndelFuncParam[2],IndelFuncParam[3]};
      //fprintf(stderr,"adding stochastic indels with parameters %f \t %f \t %f \t %f\n",pars[0],pars[1],pars[2],pars[3]);

      int ops[2] ={0,0};

      if(pars[1] == 0){
        //only potential insertions
        if(mrand_pop(mr)<pars[0]){
          if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
            add_indel_amplicon_fa(mr,&FQseq->seq,pars,ops);
          }
          else if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
            add_indel_amplicon_fqbam(mr,&FQseq->seq,&FQseq->qual,pars,ops,ErrProbTypeOffset);
          }
        }
        else{
          continue;
        }
      }
      else if(pars[0] == 0){
        //only potential deletions
        if(mrand_pop(mr)<pars[1]){
          if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
            add_indel_amplicon_fa(mr,&FQseq->seq,pars,ops);
          }
          else if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
            add_indel_amplicon_fqbam(mr,&FQseq->seq,&FQseq->qual,pars,ops,ErrProbTypeOffset);
          }
        }
        else{
          continue;
        }
      }
      else if(mrand_pop(mr)<pars[0] && mrand_pop(mr)<pars[1]){
        if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
          add_indel_amplicon_fa(mr,&FQseq->seq,pars,ops);
        }
        else if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
          add_indel_amplicon_fqbam(mr,&FQseq->seq,&FQseq->qual,pars,ops,ErrProbTypeOffset);
        }
      }

      //fprintf(stderr,"done adding insertions sequence for read \n\t%s\n\t%s\n\t%s\t\n with sizes seq %d \t qual %d\n",FQseq->name.s,FQseq->seq.s,FQseq->qual.s,FQseq->seq.l,FQseq->qual.l);
      if (ops[0] > 0 && ops[1] == 0){
        has_indels = 1;
      }
      else if (ops[0] == 0 && ops[1] > 0){
        has_indels = 2;
      }
      else if (ops[0] > 0 && ops[1] > 0){
        has_indels = 3;
      }

      //fprintf(stderr,"indel value %d\n",has_indels);
    }

    //kseq_read(seq);
    kstring_t thread_out = {0, 0, NULL};  // Initialize kstring_t for the formatted output
    
    if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
      ksprintf(&thread_out, ">%s_mod%d%d%d\n%s\n",FQseq->name.s,ReadDeam,FragMisMatch,has_indels,FQseq->seq.s);
    }
    else if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
      if(filetype == 0){
        // input amplicon file is fasta format without quality scores, thus we need to generate an artificial quality score sequence
        fprintf(stderr,"Warning without the nucleotide quality scores present in the input fasta file, an artifical quality score of 40 is created for each nucleotide\n");
        kstring_t qual;
        qual.s = NULL; qual.l=qual.m=0;
        qual.s = (char*)malloc(( FQseq->seq.l + 1) * sizeof(char));
        qual.m = FQseq->seq.m + 1;  // Set the maximum allocated length
        qual.l = FQseq->seq.l;   

        memset(qual.s, 'I', FQseq->seq.l);
        qual.s[FQseq->seq.l] = '\0'; //create the null termination in the end
        /*for (int i = 0; i < FQseq->seq.l; i++) {
          fprintf(stderr,"I VALUE %d \t %c \t %c \n",i,FQseq->seq.s[i]);
          //int quality_val = (int)(mrand_pop(mr) * 39 + 1);
          qual.s[i] = 'I' ;//(char)(quality_val + ErrProbTypeOffset);
          std::cout << qual.s[i] << std::endl;
        }*/
        ksprintf(&thread_out, "@%s_mod%d%d%d\n%s\n+\n%s\n",FQseq->name.s,ReadDeam,FragMisMatch,has_indels,FQseq->seq.s,qual.s);
        free(qual.s);
      }
      else{
        ksprintf(&thread_out, "@%s_mod%d%d%d\n%s\n+\n%s\n",FQseq->name.s,ReadDeam,FragMisMatch,has_indels,FQseq->seq.s, FQseq->qual.s);
      }
    }
    else if (amp_thread_struct->OutFormat==samT || amp_thread_struct->OutFormat==bamT){
      // consider adding ba format to simply store the sequence information

      fprintf(stderr,"Warning: NGSNGS amplicon mode on fasta or fastq files without alignment information. NGSNGS are unable to store the sequence reads in a Sequence Alignment/Map format, try using fa,fasta,fasta.gz,fa.gz,fq,fastq,fq.gz,fastq.gz output format\n");
      exit(1);
    }
    //fprintf(stderr,"TID%d_line%d_%s\t%s\t+\t%s\n", threadid, startLine,seq->name.s, seq->seq.s, seq->qual.s);
    pthread_mutex_lock(&amplicon_write_mutex);
    if (bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) < 0) {
      fprintf(stderr, "Error writing to output file in thread %d\n", threadid);
    }
    //assert(bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) != 0);
    pthread_mutex_unlock(&amplicon_write_mutex);

    free(thread_out.s);
    startLine++;
    //fprintf(stderr,"thread %d \t startline %d\n",threadid,startLine);
  }

  kseq_destroy(FQseq);

  return NULL;
}

void* ProcessBAM(void* args){
  struct_for_amplicon_threads* amp_thread_struct = (struct_for_amplicon_threads*)args;

  BGZF* amplicon_in_fp = amp_thread_struct->amplicon_in_fp;
  BGZF* amplicon_out_fp = amp_thread_struct->amplicon_out_fp;
  samFile* amplicon_in_sam = amp_thread_struct->amplicon_in_sam;
  samFile* amplicon_out_sam = amp_thread_struct->amplicon_out_sam; 
  bam_hdr_t *hdr = amp_thread_struct->hdr;

  int startLine = amp_thread_struct->startLine;
  int endLine = amp_thread_struct->endLine;
  int threadid = amp_thread_struct->threadid;

  int filetype = amp_thread_struct->filetype;
  
  //char* Briggs = amp_thread_struct->Briggs;
  //char* IndelInputParam = amp_thread_struct->Indel;

  long int Seed = amp_thread_struct->Seed;
  //fprintf(stderr,"initialize thread %d reading chunk starting from line %d to ending line %d\n",threadid,startLine,endLine);
  
  
  // Initialize briggs parameters
  float Param[4];
  if (amp_thread_struct->Briggs != NULL){
    char* BriggsParam;
    BriggsParam = strdup(amp_thread_struct->Briggs);
    Param[0] = atof(strtok(BriggsParam,"\", \t"));
    Param[1] = atof(strtok(NULL,"\", \t"));
    Param[2] = atof(strtok(NULL,"\", \t"));
    Param[3] = atof(strtok(NULL,"\", \t"));
      
    free(BriggsParam);
  }

  // initialize indel parameters
  float IndelFuncParam[4];
  if (amp_thread_struct->Indel != NULL){
    char* IndelInputParam = strdup(amp_thread_struct->Indel);
    IndelFuncParam[0] = atof(strtok(IndelInputParam,"\", \t"));
    IndelFuncParam[1] = atof(strtok(NULL,"\", \t"));
    IndelFuncParam[2] = atof(strtok(NULL,"\", \t"));
    IndelFuncParam[3] = atof(strtok(NULL,"\", \t"));
    
    free(IndelInputParam); 
  }

  //generating kstring for potential records of the stochastic indels
  char INDEL_INFO[1024];  

  // Count the number of processed reads  fprintf(stderr,"\t-> Number of reads generated by thread %d is %zu \n",struct_obj->threadno,localread);
  size_t moduloread = amp_thread_struct->moduloread;
  size_t localread = 0;
  size_t current_reads_atom = 0;

  //allocate the random generator
  mrand_t *mr = mrand_alloc(amp_thread_struct->rng_type,Seed);

  int ErrProbTypeOffset = 0;
  if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
    ErrProbTypeOffset=33; //quality score offset, 33 only for fastq it should be 0 for bam
  }

  // Initialize an alignment
  bam1_t *aln = bam_init1();

  /*size_t l_qname, const char *qname,
  uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
  size_t n_cigar, const uint32_t *cigar,
  int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
  size_t l_seq, const char *seq, const char *qual,
  size_t l_aux*/

  kstring_t Sequence;
  Sequence.s=NULL;
  Sequence.l=Sequence.m=0; // Initialize a kstring_t structure
  kstring_t Quality;
  Quality.s=NULL;
  Quality.l=Quality.m=0; // Initialize a kstring_t structure

  while (sam_read1(amplicon_in_sam, hdr, aln) >= 0) {
    size_t l_qname = aln->core.l_qname;
    const char* qname = bam_get_qname(aln);
    uint32_t flag = aln->core.flag;
    int32_t tid = aln->core.tid; 
    hts_pos_t pos = aln->core.pos;
 		uint8_t mapQ = aln->core.qual ;
    size_t n_cigar = aln->core.n_cigar;
    const uint32_t *cigar = bam_get_cigar(aln);
    int32_t mtid = aln->core.mtid;
    hts_pos_t mpos = aln->core.mpos;
    hts_pos_t isize = aln->core.isize;
    size_t l_seq = aln->core.l_qseq;
    uint8_t *seq = bam_get_seq(aln);
    uint8_t *qual = bam_get_qual(aln);
    size_t l_aux = bam_get_l_aux(aln);
    uint8_t *aux = bam_get_aux(aln);

		char *chr = hdr->target_name[aln->core.tid] ; //contig name (chromosome)
    //char *cigar_str =  PrintCigarBamSet1(n_cigar,cigar);
    //Create sequence and quality strings
    CreateSeqQualKString(aln, &Sequence, &Quality,ErrProbTypeOffset);

    //Reads aligning to reverse strand is in the orientation of the reference genome, therefore they need to be reverse complemented first
    if((aln->core.flag&BAM_FREVERSE) != 0){
      //fprintf(stderr,"qname is %s\tflag is %d\n",qname,flag);
      //fprintf(stderr,"sequence before \t %s \n",Sequence.s);
      ReversComplement_k(&Sequence);
      //fprintf(stderr,"sequence after \t %s \n",Sequence.s);
    }

    /*
    fprintf(stderr,"Sequence 1 %s\n",Sequence.s);

    SimBriggsModel_amplicon(&Sequence,0.024,0.36,0.68,0.0097,mr);
    fprintf(stderr,"Sequence 2 %s\n",Sequence.s);
    */
    localread++;
    current_reads_atom++;
    //printing out every tenth of the runtime
    if (current_reads_atom > 1 && current_reads_atom%moduloread == 0)
      fprintf(stderr,"\t-> Processed %zu reads with a current total of %zu\n",moduloread,current_reads_atom);

    // Sequence alteration integers
    int FragMisMatch = 0;
    int has_indels = 0;
    int ReadDeam=0;
    
    // deamination
    if (amp_thread_struct->Briggs != NULL){
      int strand = mrand_pop(mr)>0.5?0:1;
      ReadDeam = SimBriggsModel_amplicon(&Sequence,Param[0],Param[1],Param[2],Param[3], mr);
    }

    // Mismatch matrix input file
    if(amp_thread_struct->doMisMatchErr > 0){
      //fprintf(stderr,"INSIDE mf\n");    
      FragMisMatch = MisMatchFile_kstring(&Sequence,mr,amp_thread_struct->MisMatch,amp_thread_struct->MisLength);
      //fprintf(stderr,"FragMisMatch val %d \n",FragMisMatch);
    }

    // Stochastic structural variation model    
    if(amp_thread_struct->Indel != NULL){
      double pars[4] = {IndelFuncParam[0],IndelFuncParam[1],IndelFuncParam[2],IndelFuncParam[3]};
      //fprintf(stderr,"adding stochastic indels with parameters %f \t %f \t %f \t %f\n",pars[0],pars[1],pars[2],pars[3]);

      int ops[2] ={0,0};

      if(pars[1] == 0){
        //only potential insertions
        if(mrand_pop(mr)<pars[0]){
          if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
            add_indel_amplicon_fa(mr,&Sequence,pars,ops);
          }
          else{
            add_indel_amplicon_fqbam(mr,&Sequence,&Quality,pars,ops,ErrProbTypeOffset);
          }
        }
        else{
          continue;
        }
      }
      else if(pars[0] == 0){
        //only potential deletions
        if(mrand_pop(mr)<pars[1]){
          if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
            add_indel_amplicon_fa(mr,&Sequence,pars,ops);
          }
          else{
            add_indel_amplicon_fqbam(mr,&Sequence,&Quality,pars,ops,ErrProbTypeOffset);
          }
        }
        else{
          continue;
        }
      }
      else if(mrand_pop(mr)<pars[0] && mrand_pop(mr)<pars[1]){
        if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
          add_indel_amplicon_fa(mr,&Sequence,pars,ops);
        }
        else {
          add_indel_amplicon_fqbam(mr,&Sequence,&Quality,pars,ops,ErrProbTypeOffset);
        }
      }

      //fprintf(stderr,"done adding insertions sequence for read \n\t%s\n\t%s\n\t%s\t\n with sizes seq %d \t qual %d\n",FQseq->name.s,FQseq->seq.s,FQseq->qual.s,FQseq->seq.l,FQseq->qual.l);
      if (ops[0] > 0 && ops[1] == 0){
        has_indels = 1;
      }
      else if (ops[0] == 0 && ops[1] > 0){
        has_indels = 2;
      }
      else if (ops[0] > 0 && ops[1] > 0){
        has_indels = 3;
      }

      //fprintf(stderr,"indel value %d\n",has_indels);
    }


    if (amp_thread_struct->OutFormat==samT || amp_thread_struct->OutFormat==bamT){
      bam1_t *aln_out = bam_init1();
      //fprintf(stderr,"Sequence 2 %s\n",Sequence.s);
      
      size_t qnamelen = snprintf(NULL, 0, ">%s_mod%d%d%d", qname, ReadDeam, FragMisMatch, has_indels);
      char* formatted_qname = (char*)malloc(qnamelen + 1);
      sprintf(formatted_qname, ">%s_mod%d%d%d", qname, ReadDeam, FragMisMatch, has_indels);

      /*size_t l_qname, const char *qname,
      uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
      size_t n_cigar, const uint32_t *cigar,
      int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
      size_t l_seq, const char *seq, const char *qual,
      size_t l_aux*/

      //Reads aligning to reverse strand is in the orientation of the reference genome, before modifications they were reverse complemented, so now with bam output we change the orienation back to the reference genome
      if((aln->core.flag&BAM_FREVERSE) != 0){
        ReversComplement_k(&Sequence);
      }

      bam_set1(aln_out,qnamelen,formatted_qname,flag,tid,pos,mapQ,n_cigar,cigar,mtid,mpos,isize,Sequence.l,Sequence.s,Quality.s,l_aux);
      free(formatted_qname);

      pthread_mutex_lock(&amplicon_write_mutex);
      assert(sam_write1(amplicon_out_sam, hdr, aln_out) >=0);
      pthread_mutex_unlock(&amplicon_write_mutex);

      bam_destroy1(aln_out);
    }
    else{
      kstring_t thread_out = {0, 0, NULL};  // Initialize kstring_t for the formatted output
      
      if (amp_thread_struct->OutFormat==faT || amp_thread_struct->OutFormat==fagzT){
        ksprintf(&thread_out, ">%s_mod%d%d%d\n%s\n",qname,ReadDeam,FragMisMatch,has_indels,Sequence.s);
      }
      else if (amp_thread_struct->OutFormat==fqT || amp_thread_struct->OutFormat==fqgzT){
        ksprintf(&thread_out, "@%s_mod%d%d%d\n%s\n+\n%s\n",qname,ReadDeam,FragMisMatch,has_indels,Sequence.s, Quality.s);
      }
      //fprintf(stderr,"TID%d_line%d_%s\t%s\t+\t%s\n", threadid, startLine,seq->name.s, seq->seq.s, seq->qual.s);
      pthread_mutex_lock(&amplicon_write_mutex);
      if (bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) < 0) {
        fprintf(stderr, "Error writing to output file in thread %d\n", threadid);
      }
      //assert(bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) != 0);
      pthread_mutex_unlock(&amplicon_write_mutex);

      free(thread_out.s);
    }

    startLine++;

    /*
    
    fprintf(stderr,"queryname %s\n",bam_get_qname(aln));
    fprintf(stderr,"qualitystring %s\n",Quality.s);
    fprintf(stderr,"Sequence %s\n",Sequence.s);
    fprintf(stderr,"Print mapping quality score %d \n",mapQ);
    fprintf(stderr,"chromosome %s\n",chr);
    fprintf(stderr,"Position %d\n",pos);
    fprintf(stderr,"CIGAR %s\n",cigar_str);
    free(cigar_str);
    */
  }

  // Free allocated memory outside the loop
  free(Sequence.s);
  free(Quality.s);

  // Clean up
  bam_destroy1(aln);

  return NULL;
}


#ifdef __WITH_MAIN_AMPLICON__

int main(int argc,char **argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    AmpliconHelpPage(stderr);
    return 0;
  }
  else if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0))){
    fprintf(stderr,"\t-> ngsngs version %s: %s (htslib: %s) build(%s %s)\n",NGSNGS_RELEASE,NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    return 0;
  }
  else{
    mypars = amplicongetpars(argc,argv);
    if(mypars==NULL)
      return 1;
    
    fprintf(stderr,"\n\t-> ngsngs version %s: %s (htslib: %s) build(%s %s)\n",NGSNGS_RELEASE,NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stderr,"\t-> My commmand: %s\n",mypars->CommandRun);

    clock_t t = clock();
    time_t t2 = time(NULL);

    fprintf(stderr,"\t-> Seed is provided (-s): %ld\n",mypars->Seed/1000);

    //Initate seed
    if (mypars->rng_type == -1){
      #if defined(__linux__) || defined(__unix__)
        mypars->rng_type = 0;
      #elif defined(__APPLE__) || defined(__MACH__)
        mypars->rng_type = 3;
        //when 0 it will have problems with drand48 reentrant, will default to erand48 (MacroRandType 3)
      #else
      #   error "Unknown compiler"
      #endif
    }

    const char* Amplicon_in = mypars->Amplicon_in_pars;
    const char* dot = strrchr(Amplicon_in, '.');
    // Extract the file extension
    const char* extension = dot + 1;
    int filetype;
    // Print the file extension for debugging purposes
        
    if (strcmp(extension, "gz") == 0) {
      BGZF* fp_tmp = bgzf_open(Amplicon_in, "r");

      char ch;
      int result = bgzf_read(fp_tmp, &ch, 1);
      if (ch == '>'){
        filetype = 0;
        fprintf(stderr,"\t-> Amplicon mode is chosen with a compressed fasta input file\n");
      }
      else if (ch == '@'){
        fprintf(stderr,"\t-> Amplicon mode is chosen with a compressed fastq input file\n");
        filetype = 1;
      }
      else {
        fprintf(stderr,"Unknown file format");
        bgzf_close(fp_tmp);
        exit(1);
      }
      bgzf_close(fp_tmp);
    }
    else if (strcmp(extension, "fa") == 0 || strcmp(extension, "fas") == 0 || strcmp(extension, "fasta") == 0){
      filetype = 0;
      fprintf(stderr,"\t-> Amplicon mode is chosen with a fasta input file\n");
    }
    else if (strcmp(extension, "fq") == 0 || strcmp(extension, "fastq") == 0){
      fprintf(stderr,"\t-> Amplicon mode is chosen with a fastq input file\n");
      filetype = 1;
    }
    else if (strcmp(extension, "sam") == 0 || strcmp(extension, "bam") == 0){
      fprintf(stderr,"\t-> Amplicon mode is chosen with a Sequence Alignment/Map format input file\n");
      filetype = 2;
    }
    else {
      fprintf(stderr,"Unknown file format");
      exit(1);
    }

    if(mypars->OutFormat == faT){fprintf(stderr,"\t-> Amplicon mode is chosen with a fasta output file\n");}
    else if(mypars->OutFormat == fagzT){fprintf(stderr,"\t-> Amplicon mode is chosen with a compressed fasta output file\n");}
    else if(mypars->OutFormat == fqT){fprintf(stderr,"\t-> Amplicon mode is chosen with a compressed fastq output file\n");}
    else if(mypars->OutFormat == fqgzT){fprintf(stderr,"\t-> Amplicon mode is chosen with a compressed fastq output file\n");}
    else if(mypars->OutFormat == samT || mypars->OutFormat == bamT){fprintf(stderr,"\t-> Amplicon mode is chosen with a Sequence Alignment/Map format output file\n");}

    if (filetype < 2){
      // fasta or fastq

      // Count the total number of lines in the file
      size_t totalLines = 0;

      BGZF* fp_tmp = bgzf_open(mypars->Amplicon_in_pars, "r");

      kstring_t linecounttmp;linecounttmp.s=NULL;linecounttmp.l=linecounttmp.m=0; // Initialize a kstring_t structure
      while (bgzf_getline(fp_tmp, '\n', &linecounttmp) != -1){totalLines++;}
      free(linecounttmp.s);
      bgzf_close(fp_tmp);
      
      int readlinestructure;
      if(filetype < 1){//fasta
        readlinestructure = 2;
      }
      else{ //fastq
        readlinestructure = 4;
      }
      // Calculate the number of lines per thread
      size_t no_reads = totalLines/readlinestructure;
      size_t linesPerThread = totalLines / mypars->Threads;

      int modulovalue;
      if (no_reads > 1000000){
        modulovalue = 10;
      }
      else{ 
        modulovalue = 1;
      }
      size_t moduloread = no_reads/modulovalue;

      //fprintf(stderr,"Total lines %zu \t Total reads %zu \t lines pr threads %zu \n",totalLines,no_reads,mypars->Threads);  

      AmpliconThreadInitialize(mypars->OutFormat,mypars->Amplicon_in_pars,filetype,mypars->Amplicon_out_pars,mypars->SubProfile,mypars->Threads,mypars->BriggsBiotin,mypars->Indel,mypars->Seed,mypars->rng_type,moduloread,totalLines,linesPerThread);
    }
    else if (filetype == 2){
      // Count the total number of lines in the file
      size_t totalLines = 0;
      
      samFile *AmpliconBam = hts_open(mypars->Amplicon_in_pars, "r");
      bam_hdr_t *hdr = sam_hdr_read(AmpliconBam);  // Read the header
      bam1_t *aln = bam_init1();   // Initialize an alignment
      while(sam_read1(AmpliconBam, hdr, aln) >= 0){totalLines++;}
      // Clean up
      bam_destroy1(aln);
      sam_hdr_destroy(hdr);
      sam_close(AmpliconBam);
      //fprintf(stderr,"total lines in bam file %d\n",totalLines);
      size_t linesPerThread = totalLines / mypars->Threads;

      int modulovalue;
      if (totalLines > 100){
        modulovalue = 10;
      }
      else{ 
        modulovalue = 1;
      }
      size_t moduloread = totalLines/modulovalue;

      //insert thread initialization for bam with different read in function
      AmpliconThreadInitialize(mypars->OutFormat,mypars->Amplicon_in_pars,filetype,mypars->Amplicon_out_pars,mypars->SubProfile,mypars->Threads,mypars->BriggsBiotin,mypars->Indel,mypars->Seed,mypars->rng_type,moduloread,totalLines,linesPerThread);
      /*ProcessBAM(mypars->Amplicon_in_pars);
      exit(1);*/
    }

    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  amplicongetpars_destroy(mypars);
  return 0;
}
#endif

/*
g++ Amplicon.cpp -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto ../mrand.o ../Briggs.o ../NtSubModels.o ../add_indels.o ../../htslib/libhts.a -D __WITH_MAIN__ -o Amplicon

./Amplicon --amplicon Amplicon_in.fq -m b,0.024,0.36,0.68,0.0097 -mf ../Test_Examples/MisincorpFile.txt --output Amplicon_out.fq

Deletions
./Amplicon --amplicon Amplicon_in.fq -indel 0.0,0.5,0.0,0.9 --output Amplicon_out.fq -t 1

Insertions
./Amplicon --amplicon Amplicon_in.fq -indel 0.05,0.0,0.1,0.0 --output Amplicon_out.fq -t 1

./Amplicon --amplicon Amplicon_in.fq --output Amplicon_out.fq --threads 2


*/
