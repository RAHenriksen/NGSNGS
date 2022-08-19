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

  mypars->CommandRun = (char*) calloc(1024,1);
  char *Command = mypars->CommandRun;
  const char *first = "./ngsngs ";
  strcpy(Command,first);

  const char* readstr;
  const char* bufstr;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Reference = strdup(*(++argv));
      strcat(Command,mypars->Reference); strcat(Command," ");
    }
    else if(strcasecmp("-bcf",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->vcffile = strdup(*(++argv));
      strcat(Command,mypars->vcffile); strcat(Command," ");
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->SamplThreads = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->SamplThreads < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->CompressThreads = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->CompressThreads < 0){ErrMsg(9.0);}
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      readstr = strdup(*(++argv));
      sscanf(readstr, "%lu",&mypars->nreads);
      strcat(Command,*argv); strcat(Command," ");
      free((char*)readstr);
      if (mypars->nreads <= 0){ErrMsg(2.4);}
    }
    else if(strcasecmp("-bl",*argv)==0 || strcasecmp("--bufferlength",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      bufstr = strdup(*(++argv));
      sscanf(bufstr, "%lu",&mypars->KstrBuf);
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-c",*argv)==0 || strcasecmp("--cov",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->coverage = atof(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->coverage < 0.0){ErrMsg(2.2);}
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->OutName = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Glob_seed = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      char * tok = *(++argv);
      strcat(Command,*argv); strcat(Command," ");
      
      if(strcasecmp("SE",tok)==0)
	      mypars->seq_type = SE;
      else if(strcasecmp("PE",tok)==0)
	      mypars->seq_type = PE;
      if(mypars->seq_type==unknownTT)
	      ErrMsg(6.5);
    }
    else if(strcasecmp("-a1",*argv)==0 || strcasecmp("--adapter1",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Adapter1 = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-a2",*argv)==0 || strcasecmp("--adapter2",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Adapter2 = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-q1",*argv)==0 || strcasecmp("--quality1",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->QualProfile1 = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-q2",*argv)==0 || strcasecmp("--quality2",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->QualProfile2 = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-mf",*argv)==0 || strcasecmp("--mismatch",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->SubProfile = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-ne",*argv)==0 || strcasecmp("--noerror",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->DoSeqErr = 0;//"F";
    }
    else if(strcasecmp("-na",*argv)==0 || strcasecmp("--noalign",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->NoAlign = 1;
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      ++argv;
      strcat(Command,*argv); strcat(Command," ");
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
      strcat(Command,*argv); strcat(Command," ");
      mypars->Briggs = strdup(*(++argv)); //double nv, double lambda, double delta_s, double delta -> 0.024,0.36,0.68,0.0097
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Length = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->Length < 0.0){ErrMsg(3.2);}
    }
    else if(strcasecmp("-lf",*argv)==0 || strcasecmp("--lengthfile",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->LengthFile = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-ld",*argv)==0 || strcasecmp("--lengthdist",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->LengthDist = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-chr",*argv)==0 || strcasecmp("--chromosomes",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Chromosomes = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-p",*argv)==0 || strcasecmp("--poly",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Poly = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if(strcasecmp("A",mypars->Poly)!=0 && 
      strcasecmp("G",mypars->Poly)!=0 &&
      strcasecmp("C",mypars->Poly)!=0 &&
      strcasecmp("T",mypars->Poly)!=0 &&
      strcasecmp("N",mypars->Poly)!=0){ErrMsg(10.0);}
    }
    else if(strcasecmp("-id",*argv)==0 || strcasecmp("--indiv",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->HeaderIndiv = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-rng",*argv)==0 || strcasecmp("--rand",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->rng_type = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else{
      fprintf(stderr,"Unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      return NULL;
    }
    
    ++argv;
  }
  //free((char*)bufstr);
  //free(mypars->CommandRun)
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
