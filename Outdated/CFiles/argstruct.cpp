#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>//for printing time

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <cstdlib>
#include <ctime>

#include <cstdio>
#include <cassert>
#include <cstdint>

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>        
#include <atomic>
#include <vector>

#include <getopt.h>
#include "getopt_func.h"

typedef struct{
  int reads;
  int threads;
  char *out;
  const char *fastafile;

  const char* Adapt_flag; //optional
  int seeds;
}argStruct;

int print_pars(argStruct *p,FILE *fp){
  fprintf(fp,"Parsed arguments \t out: %s Reads: %d\n",p->out,p->reads);
  return 0;
}

int HelpPage(FILE *fp){
  fprintf(fp,"./ngsngs [-type1,type2] [-typeX int,-out prefix]\n");
  return 0;
}

//returtype navn (par1,par2,par3)

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->reads = -1;
  mypars->out = strdup("output");

  while(*argv){
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->fastafile = strdup(*(++argv));
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->reads = atoi(*(++argv));
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->threads = atoi(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--out",*argv)==0){
      mypars->out = strdup(*(++argv));
    }
    // optional arguments 
    else if(strcasecmp("-a",*argv)==0 || strcasecmp("--adapter",*argv)==0){
      mypars->Adapt_flag = strdup(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->seeds = atoi(*(++argv));
    }
    else{
      HelpPage(stderr);
      ++argv;
    } //hvordan tager vi højde for argument som ikke findes?
    ++argv;
  }
  return mypars;
}

int main(int argc,char**argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                        strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    fprintf(stderr,"help page statement\n");
    HelpPage(stderr);
    return 0;
  }
  else{
    mypars = getpars(argc,argv);
  }
  
  Filecreate(mypars->out,mypars->reads);  
  //print_type(mypars->reads);
  print_pars(mypars,stdout);
  return 0;
}