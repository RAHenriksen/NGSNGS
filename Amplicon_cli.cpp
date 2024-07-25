#include <htslib/kstring.h>
#include "Amplicon_cli.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <cassert>

argStruct *amplicongetpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;

  mypars->Amplicon_in_pars = NULL;
  mypars->Amplicon_out_pars = NULL;
  mypars->BriggsBiotin = NULL;
  mypars->Threads = 1;
  mypars->Seed = (int) time(NULL);
  mypars->rng_type = -1;
  // 3) misincorporation matrix
  mypars->SubProfile = NULL;
  mypars->Indel = NULL;

  //quality
  mypars->fixqual = 0;
  mypars->DoSeqErr = 1;
  mypars->QualProfile = NULL;
  
  mypars->OutFormat = unknownT;
  
  mypars->CommandRun = NULL;
  kstring_t kstr;kstr.s=NULL;kstr.l=kstr.m=0;
  for(int i=0;i<argc;i++)
    ksprintf(&kstr,"%s ",argv[i]);
  mypars->CommandRun = kstr.s;

  ++argv;
  while(*argv){
    if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      ++argv;
      char *tok = *argv;
      if(strcasecmp("fa",tok)==0 || strcasecmp("fasta",tok)==0)
	      mypars->OutFormat = faT;
      if(strcasecmp("fa.gz",tok)==0 || strcasecmp("fasta.gz",tok)==0)
	      mypars->OutFormat = fagzT;
      if(strcasecmp("fq",tok)==0 || strcasecmp("fastq",tok)==0)
	      mypars->OutFormat = fqT;
      if(strcasecmp("fq.gz",tok)==0 || strcasecmp("fastq.gz",tok)==0)
	      mypars->OutFormat = fqgzT;
      if(strcasecmp("sam",tok)==0)
	      mypars->OutFormat = samT;
      if(strcasecmp("bam",tok)==0)
	      mypars->OutFormat = bamT;
      if(mypars->OutFormat==unknownT){
	      fprintf(stderr,"\nNext Generation Simulator for Next Generator Sequencing Data Amplicon\nWarning:\n");
        ErrMsg(7.0);
      }
    }
    else if(strcasecmp("-a",*argv)==0 || strcasecmp("--amplicon",*argv)==0){
      mypars->Amplicon_in_pars = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->Amplicon_out_pars = strdup(*(++argv));
    }
    else if(strcasecmp("-m",*argv)==0 || strcasecmp("--model",*argv)==0){
      ++argv;
      char *tok = *argv;
      char* ModelString = strdup(tok);
      char* BriggsModel;
      BriggsModel = strtok(ModelString,",");
      char* ModelParam =  strdup(strtok (NULL, ""));
      if(strcasecmp("b",BriggsModel)==0 || strcasecmp("briggs",BriggsModel)==0 || strcasecmp("b7",BriggsModel)==0 || strcasecmp("briggs07",BriggsModel)==0)
	      mypars->BriggsBiotin = ModelParam;

      free(ModelString);
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->Threads = atoi(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Seed = atol(*(++argv))*1000;
    }
    else if(strcasecmp("-mf",*argv)==0 || strcasecmp("--mismatch",*argv)==0){
      mypars->SubProfile = strdup(*(++argv));
    }
    else if(strcasecmp("-indel",*argv)==0){
      mypars->Indel = strdup(*(++argv));
    }
    else if(strcasecmp("-rng",*argv)==0 || strcasecmp("--rand",*argv)==0){
      mypars->rng_type = atoi(*(++argv));
    }
    else if(strcasecmp("-qs",*argv)==0 || strcasecmp("--qualityscore",*argv)==0){
      mypars->fixqual = atoi(*(++argv));
    }
    else if(strcasecmp("-ne",*argv)==0 || strcasecmp("--noerror",*argv)==0){
      mypars->DoSeqErr = 0;//"F";
    }
    else if(strcasecmp("-q",*argv)==0 || strcasecmp("--quality",*argv)==0){
      mypars->QualProfile = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

void amplicongetpars_destroy(argStruct *mypars){
  free(mypars->CommandRun);
  free(mypars->BriggsBiotin);
  free(mypars->SubProfile);
  free(mypars->Indel);
  free(mypars->QualProfile);
  delete mypars;
}
