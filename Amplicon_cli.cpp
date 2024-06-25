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

  ++argv;
  while(*argv){
    if(strcasecmp("-a",*argv)==0 || strcasecmp("--amplicon",*argv)==0){
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
      mypars->Seed = atol(*(++argv));
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
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}
