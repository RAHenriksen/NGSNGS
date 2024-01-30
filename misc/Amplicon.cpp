#include "../mrand.h"
#include "../fasta_sampler.h"
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

#define MAX_LINE_LENGTH 2048

pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

// Define parameters

typedef struct{
  const char* Amplicon_in_pars;
  const char* Amplicon_out_pars;
  int Threads;
}argStruct;

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;

  mypars->Amplicon_in_pars = NULL;
  mypars->Amplicon_out_pars = NULL;
  mypars->Threads = 1;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-a",*argv)==0 || strcasecmp("--amplicon",*argv)==0){
      mypars->Amplicon_in_pars = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->Amplicon_out_pars = strdup(*(++argv));
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->Threads = atoi(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

// Define a struct to store the information from a FastQ record
typedef struct {
    char* header;
    char* sequence;
    char* quality;
} FastQRecord;

typedef struct {
  BGZF* amplicon_in_fp;
  int startLine;
  int endLine;
  int threadid;
  BGZF* amplicon_out_fp;  // Include the output file in the thread struct
} struct_for_amplicon_threads;


KSEQ_INIT(BGZF*, bgzf_read)

void* ProcessChunk(void* args) {
  struct_for_amplicon_threads* amp_thread_struct = (struct_for_amplicon_threads*)args;
  BGZF* amplicon_in_fp = amp_thread_struct->amplicon_in_fp;
  BGZF* amplicon_out_fp = amp_thread_struct->amplicon_out_fp;
  int startLine = amp_thread_struct->startLine;
  int endLine = amp_thread_struct->endLine;
  int threadid = amp_thread_struct->threadid;
  
  fprintf(stderr,"initialize thread %d reading chunk starting from line %d to ending line %d\n",threadid,startLine,endLine);
  // Seek to the appropriate starting line
  kseq_t *seq = kseq_init(amplicon_in_fp);
  while (startLine <= endLine) {
    if (kseq_read(seq) < 0) {
      fprintf(stderr, "Error reading sequence in thread %d\n", threadid);
      break;    
    }
    //kseq_read(seq);
    //fprintf(stderr,"thread %d for current line %d and end line %d\n",threadid,currentLine,endLine);
    kstring_t thread_out = {0, 0, NULL};  // Initialize kstring_t for the formatted output
    ksprintf(&thread_out, "TID%d_line%d_%s\n%s\n+\n%s\n", threadid, startLine,seq->name.s, seq->seq.s, seq->qual.s);
    fprintf(stderr,"TID%d_line%d_%s\t%s\t+\t%s\n", threadid, startLine,seq->name.s, seq->seq.s, seq->qual.s);
    pthread_mutex_lock(&write_mutex);
    if (bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) < 0) {
      fprintf(stderr, "Error writing to output file in thread %d\n", threadid);
    }
    //assert(bgzf_write(amplicon_out_fp, thread_out.s, thread_out.l) != 0);
    pthread_mutex_unlock(&write_mutex);

    free(thread_out.s);
    startLine++;
    fprintf(stderr,"thread %d \t startline %d\n",threadid,startLine);
  }

  kseq_destroy(seq);
  return NULL;
}

#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  argStruct *mypars = NULL;
  mypars = getpars(argc,argv);
  
  const char* Amplicon_in_fp = mypars->Amplicon_in_pars;
  const char* Amplicon_out_fp = mypars->Amplicon_out_pars;
  int threads = mypars->Threads;

  //const char* fastqFile = "Amplicon_in.fq";
  
  // Count the total number of lines in the file
  long long totalLines = 0;
  BGZF* fp_tmp = bgzf_open(Amplicon_in_fp, "r");
  kstring_t linecounttmp;linecounttmp.s=NULL;linecounttmp.l=linecounttmp.m=0; // Initialize a kstring_t structure
  while (bgzf_getline(fp_tmp, '\n', &linecounttmp) != -1){totalLines++;}
  free(linecounttmp.s);
  // Calculate the number of lines per thread
  size_t no_reads = totalLines/4;
  size_t linesPerThread = totalLines / threads;
  bgzf_close(fp_tmp);

  printf("Total lines %zu \t Total reads %zu \t lines pr threads %zu \n",totalLines,no_reads,linesPerThread);
  
  BGZF* fp = bgzf_open(Amplicon_in_fp, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error opening the FastQ file.\n");
    return 1;
  }

  BGZF* amplicon_out = bgzf_open(Amplicon_out_fp, "wu");
  if (amplicon_out == NULL) {
    fprintf(stderr, "Error opening output file.\n");
    return 1;
  }

  // Create an array to hold thread IDs dynamically
  pthread_t* mythreads = new pthread_t[threads];

  // Create an array to hold thread arguments dynamically
  struct_for_amplicon_threads* amp_thread_struct = new struct_for_amplicon_threads[threads];

  // Create threads dynamically
  for (int i = 0; i < threads; ++i) {
    int startLine = i * linesPerThread+1;
    int endLine = (i == threads - 1) ? totalLines : startLine + linesPerThread - 1;

    fprintf(stderr, "Thread %d reads in line start %d to line end %d\n", i, startLine, endLine);

    amp_thread_struct[i].amplicon_in_fp = fp;
    amp_thread_struct[i].amplicon_out_fp = amplicon_out;
    amp_thread_struct[i].startLine = startLine;
    amp_thread_struct[i].endLine = endLine;
    amp_thread_struct[i].threadid = i;

    // Create each thread
    pthread_create(&mythreads[i], NULL, ProcessChunk, &amp_thread_struct[i]);
  }
  
  // Wait for all threads to finish
  for (int i = 0; i < threads; ++i) {
    pthread_join(mythreads[i], NULL);
    fprintf(stderr, "Thread %d finished\n", i);
  }

  // Delete allocated memory
  delete[] mythreads;
  delete[] amp_thread_struct;

  // Close the files
  bgzf_close(amplicon_out);
  bgzf_close(fp);

  fprintf(stderr, "done closing the file\n");
  
  //fprintf(stderr, "done closing the file\n");
  return 0;
}
#endif

/*
g++ Amplicon.cpp -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto ../../htslib/libhts.a -D __WITH_MAIN__ -o Amplicon

./Amplicon --amplicon Amplicon_in.fq --output Amplicon_out.fq --threads 2
*/

/*
  // CHUNK
  // Calculate the virtual offset for the beginning of line 100
    int targetLine = 100;
    int offset = bgzf_tell(fp);  // Get the current offset
    int line = 0;

    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0 && line < targetLine) {
        offset = bgzf_tell(fp);
        line++;
    }

    // Use bgzf_seek to set the file position to the beginning of line 100
    bgzf_seek(fp, offset, SEEK_SET);

    // Read and process lines 100 to 200
    while (kseq_read(seq) >= 0 && line < targetLine + 100) {
        printf("Header: %s\n", seq->name.s);
        printf("Sequence: %s\n", seq->seq.s);
        printf("Quality: %s\n", seq->qual.s);

        line++;
    }

  kseq_destroy(seq);
  bgzf_close(fp);
  */



  /*
  READ IN
  kseq_t *seq = kseq_init(fp);
  while (kseq_read(seq) >= 0) {
    FastQRecord record;
    record.header = strdup(seq->name.s);
    record.sequence = strdup(seq->seq.s);
    record.quality = strdup(seq->qual.s);

    // Process the current record as needed
    printf("Header: %s\n", record.header);
    printf("Sequence: %s\n", record.sequence);
    printf("Quality: %s\n", record.quality);

    // Free memory allocated by strdup
    free(record.header);
    free(record.sequence);
    free(record.quality);
    printf("done freeing\n");
  }
  kseq_destroy(seq);
  bgzf_close(fp);
  */