#include <cstdio>
#include <cassert>
#include <htslib/faidx.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <typeinfo>

// I would like to create a function with TK's code since its optimal in case we wish to 
//simulate a given number of fragments

void random_seq(faidx_t *seq_ref){
  // choose a random sequence -> still ned to change it so it saves the output to a single file.
  int readlength=35;
  int nreads = 8;
  for(int i=0;i<nreads;i++){
    char buf[96];//assume maxlength for readid is 96bytes
    int whichref = lrand48() % faidx_nseq(seq_ref);
    fprintf(stderr,"\t-> Whichref: %d\n",whichref);
    const char *name = faidx_iseq(seq_ref,whichref);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"\t-> name: \'%s\' name_len: %d\n",name,name_len);

    int start = lrand48() % name_len;
    int stop = start+readlength;
    if(stop>name_len)
      stop = name_len;
    snprintf(buf,96,"%s:%d-%d",name,start,stop);
    fprintf(stderr,"buf: %s\n",buf);
    fprintf(stdout,"%s\n+\n",buf);
    char *data = fai_fetch(seq_ref,name,&name_len);

    for(int i=start;i<stop;i++)
      fprintf(stdout,"%c",data[i]);
    fprintf(stdout,"\n");

    for(int i=start;i<stop;i++)
      fprintf(stdout,"F");
    fprintf(stdout,"\n");
  }
}

int main(int argc,char **argv){

  const char *fastafile = "chr22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  
  // choosing random sequences using -> random_seq(ref);

  // is lrand48() in order to pick a random sequence if containing more?
  int whichref = lrand48() % faidx_nseq(ref);
  const char *name = faidx_iseq(ref,whichref);
  int name_len =  faidx_seq_len(ref,name);
  std::cout << name << std::endl;
  std::cout << name_len << std::endl;
  
  char *sequence = faidx_fetch_seq(ref,name,30000000,30000100,&name_len);
  std::cout << sequence << std::endl;

  return 0;
}
