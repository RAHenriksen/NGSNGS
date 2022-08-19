#include <htslib/kstring.h>
#include "NGSNGS_cli.h"

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->SamplThreads = 1;
  mypars->CompressThreads = 1;

  // generating strings for which the simulated reads will be contained
  mypars->nreads = 0;
  mypars->coverage = 0.0;
  mypars->KstrBuf=30000000;

  // The output format, output files, and structural elements for SAM outputs
  mypars->OutFormat = unknownT;
  mypars->OutName = NULL; //"output";
  mypars->HeaderIndiv=NULL;
  mypars->NoAlign=0;

  // Thread generation and sampling specific information
  mypars->Chromosomes = NULL;
  mypars->Reference = NULL;
  mypars->seq_type = unknownTT;
  mypars->Glob_seed = (int) time(NULL);
  mypars->rng_type = -1;

  // Sequence alteration models
  // 1) nucleotide quality score and sequencing errors,  
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  mypars->DoSeqErr = 1;
  // 2) briggs model
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  // 3) misincorporation matrix
  mypars->SubProfile = NULL;
  // 4) Bcf file and variation incorporation
  mypars->vcffile = NULL;

  // Fragment lengths 
  mypars->Length = 0;
  mypars->LengthFile = NULL;
  mypars->LengthDist = NULL;

  // Additional information for sequence reads
  mypars->Adapter1 = NULL;
  mypars->Adapter2 = NULL;
  mypars->Poly = NULL;

  mypars->CommandRun = NULL;
  kstring_t kstr;kstr.s=NULL;kstr.l=kstr.m=0;
  for(int i=0;i<argc;i++)
    ksprintf(&kstr,"%s ",argv[i]);
  mypars->CommandRun = kstr.s;//notice that this is allocated here but deallocated somewhere else

  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-bcf",*argv)==0){
      mypars->vcffile = strdup(*(++argv));
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->SamplThreads = atoi(*(++argv));
      if (mypars->SamplThreads < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      mypars->CompressThreads = atoi(*(++argv));
      if (mypars->CompressThreads < 0){ErrMsg(9.0);}
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->nreads = atol(*(++argv));
      if (mypars->nreads <= 0){ErrMsg(2.4);}
    }
    else if(strcasecmp("-bl",*argv)==0 || strcasecmp("--bufferlength",*argv)==0){
      mypars->KstrBuf = atol(*(++argv));
    }
    else if(strcasecmp("-c",*argv)==0 || strcasecmp("--cov",*argv)==0){
      mypars->coverage = atof(*(++argv));
      if (mypars->coverage < 0.0){ErrMsg(2.2);}
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->OutName = strdup(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv));
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      char * tok = *(++argv);
      
      if(strcasecmp("SE",tok)==0)
	      mypars->seq_type = SE;
      else if(strcasecmp("PE",tok)==0)
	      mypars->seq_type = PE;
      if(mypars->seq_type==unknownTT)
	      ErrMsg(6.5);
    }
    else if(strcasecmp("-a1",*argv)==0 || strcasecmp("--adapter1",*argv)==0){
      mypars->Adapter1 = strdup(*(++argv));
    }
    else if(strcasecmp("-a2",*argv)==0 || strcasecmp("--adapter2",*argv)==0){
      mypars->Adapter2 = strdup(*(++argv));
    }
    else if(strcasecmp("-q1",*argv)==0 || strcasecmp("--quality1",*argv)==0){
      mypars->QualProfile1 = strdup(*(++argv));
    }
    else if(strcasecmp("-q2",*argv)==0 || strcasecmp("--quality2",*argv)==0){
      mypars->QualProfile2 = strdup(*(++argv));
    }
    else if(strcasecmp("-mf",*argv)==0 || strcasecmp("--mismatch",*argv)==0){
      mypars->SubProfile = strdup(*(++argv));
    }
    else if(strcasecmp("-ne",*argv)==0 || strcasecmp("--noerror",*argv)==0){
      mypars->DoSeqErr = 0;//"F";
    }
    else if(strcasecmp("-na",*argv)==0 || strcasecmp("--noalign",*argv)==0){
      mypars->NoAlign = 1;
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      ++argv;
      char *tok = *argv;
      if(strcasecmp("fa",tok)==0)
	mypars->OutFormat = faT;
      if(strcasecmp("fa.gz",tok)==0)
	mypars->OutFormat = fagzT;
      if(strcasecmp("fq",tok)==0)
	mypars->OutFormat = fqT;
      if(strcasecmp("fq.gz",tok)==0)
	mypars->OutFormat = fqgzT;
      if(strcasecmp("sam",tok)==0)
	mypars->OutFormat = samT;
      if(strcasecmp("bam",tok)==0)
	mypars->OutFormat = bamT;
      if(strcasecmp("cram",tok)==0)
	mypars->OutFormat = cramT;
      if(mypars->OutFormat==unknownT)
	ErrMsg(7.0);
    }
    else if(strcasecmp("-b",*argv)==0 || strcasecmp("--briggs",*argv)==0){
      mypars->Briggs = strdup(*(++argv)); //double nv, double lambda, double delta_s, double delta -> 0.024,0.36,0.68,0.0097
    }
    else if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      mypars->Length = atoi(*(++argv));
      if (mypars->Length < 0.0){ErrMsg(3.2);}
    }
    else if(strcasecmp("-lf",*argv)==0 || strcasecmp("--lengthfile",*argv)==0){
      mypars->LengthFile = strdup(*(++argv));
    }
    else if(strcasecmp("-ld",*argv)==0 || strcasecmp("--lengthdist",*argv)==0){
      mypars->LengthDist = strdup(*(++argv));
    }
    else if(strcasecmp("-chr",*argv)==0 || strcasecmp("--chromosomes",*argv)==0){
      mypars->Chromosomes = strdup(*(++argv));
    }
    else if(strcasecmp("-p",*argv)==0 || strcasecmp("--poly",*argv)==0){
      mypars->Poly = strdup(*(++argv));
      if(strcasecmp("A",mypars->Poly)!=0 && 
      strcasecmp("G",mypars->Poly)!=0 &&
      strcasecmp("C",mypars->Poly)!=0 &&
      strcasecmp("T",mypars->Poly)!=0 &&
      strcasecmp("N",mypars->Poly)!=0){ErrMsg(10.0);}
    }
    else if(strcasecmp("-id",*argv)==0 || strcasecmp("--indiv",*argv)==0){
      mypars->HeaderIndiv = strdup(*(++argv));
    }
    else if(strcasecmp("-rng",*argv)==0 || strcasecmp("--rand",*argv)==0){
      mypars->rng_type = atoi(*(++argv));
    }
    else{
      fprintf(stderr,"Unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      return NULL;
    }
    
    ++argv;
  }
  return mypars;
}


void argStruct_destroy(argStruct *mypars){
  free(mypars->Reference); //-i
  free(mypars->OutName);
  free(mypars->LengthFile);
  free(mypars->LengthDist);
  
  free(mypars->CommandRun);
  if(mypars->Chromosomes)
    free(mypars->Chromosomes);

  // OPTIONAL DEALLOCATIONS
  free(mypars->vcffile);
  free(mypars->HeaderIndiv);
  free(mypars->Adapter1);
  free(mypars->Adapter2);
  free(mypars->QualProfile1);
  free(mypars->QualProfile2);
  free(mypars->SubProfile);
  free(mypars->Briggs);
  free(mypars->Poly);
  delete mypars;
}
