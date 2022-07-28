#ifndef NGSNGS_H
#define NGSNGS_H


typedef struct{
  int SamplThreads; //sampling threads, used internally by pthread_create
  int CompressThreads; //compression threads, used external by bgzf_mt and set_hts_options
  size_t nreads; // Thorfinn made mistake
  double coverage; // Coverage of breadth, estimated from rlen, flen, genomsize, superfancy
  int Glob_seed; // local seeds are computed determistly from global. Only one seed needs to be supplied
  const char *OutFormat; //fq, fq.gz, fa, fa.gz, sam, bam, cram, ultrafancy
  const char *OutName; //prefix for output, final will be determined by [OutName.OutFormat]
  const char *Reference; //full filename for reference fasta
  const char *Seq; // singleend or paired end.
  const char *Adapter1; //actual adapterseq, R1, not flipped, reversed or completemented
  const char *Adapter2; //actual adapterseq, R2, not flipped, reversed or completemented
  const char *QualProfile1;// filename for quality prof, for R1 or SE
  const char *QualProfile2;// filename for quality prof, for R2 used only in PE
  const char *SubProfile;// filename for misincorperation, typespecific and position specific
  const char *ErrorFlag;// should we not add errors?, {T,F}
  const char *Briggs; // the four briggs parameters
  int Length; // fragment length when fixed
  const char *LengthFile; //filename for distribution of frag lengths
  const char *LengthDist;//name of pdf used for simulating fraglengths (incl parameters)
  const char *Poly; // should poly be added possible values -p G or -p A
  const char *Chromosomes; //for subsetting chromosome of interested, -chr chrX
  int rand_val; //RNG type, drand48 or rand or drand48_r etc
  const char *Variant; // filename for bcf
  const char *Variant_type; //-v SNP -v INDEL
  char *VariantFlag; //bcf
  char *CommandRun; // actual command run in same order
  const char* HeaderIndiv; //samplename from VCF/BCF file
  const char* NoAlign;// This option is cool, but explaining it takes up to much space in comment
  size_t KstrBuf; // Buffer size for kstring length
}argStruct;

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->SamplThreads = 1;
  mypars->CompressThreads = 1;

  // generating strings for which the simulated reads will be contained
  mypars->nreads = 0;
  mypars->coverage = 0.0;
  mypars->KstrBuf=30000000;

  // The output format, output files, and structural elements for SAM outputs
  mypars->OutFormat = NULL; //"fa";
  mypars->OutName = NULL; //"output";
  mypars->HeaderIndiv=NULL;
  mypars->NoAlign=NULL;

  // Thread generation and sampling specific information
  mypars->Chromosomes = NULL;
  mypars->Reference = NULL;
  mypars->Seq = NULL; // "SE";
  mypars->Glob_seed = (int) time(NULL);
  mypars->rand_val = -1;

  // Sequence alteration models
  // 1) nucleotide quality score and sequencing errors,  
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  mypars->ErrorFlag = NULL;
  // 2) briggs model
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  // 3) misincorporation matrix
  mypars->SubProfile = NULL;
  // 4) Bcf file and variation incorporation
  mypars->Variant = NULL;
  mypars->Variant_type = NULL;

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
      mypars->Variant = strdup(*(++argv));
      strcat(Command,mypars->Variant); strcat(Command," ");
    }
    else if(strcasecmp("-v",*argv)==0 || strcasecmp("--variant",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->Variant_type = strdup(*(++argv));
      strcat(Command,mypars->Variant_type); strcat(Command," ");
      if(mypars->Variant == NULL){ErrMsg(13.0);}
      else if (mypars->Variant_type && strcasecmp("snp",mypars->Variant_type)!=0 && 
      strcasecmp("indel",mypars->Variant_type)!=0 &&
      strcasecmp("all",mypars->Variant_type)!=0){ErrMsg(13.5); exit(0);}      
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
      mypars->Seq = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if(strcasecmp("SE",mypars->Seq)!=0 && strcasecmp("PE",mypars->Seq)!=0){ErrMsg(6.5);} 
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
      mypars->ErrorFlag = "F";
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-na",*argv)==0 || strcasecmp("--noalign",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->NoAlign = "T";
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->OutFormat = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if(strcasecmp("fa",mypars->OutFormat)!=0 && 
      strcasecmp("fa.gz",mypars->OutFormat)!=0 &&
      strcasecmp("fq",mypars->OutFormat)!=0 &&
      strcasecmp("fq.gz",mypars->OutFormat)!=0 &&
      strcasecmp("sam",mypars->OutFormat)!=0 &&
      strcasecmp("bam",mypars->OutFormat)!=0 &&
      strcasecmp("cram",mypars->OutFormat)!=0){ErrMsg(7.0);}      
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
      mypars->rand_val = atoi(*(++argv));
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

void handler(int s);

void catchkill();

#endif

// -ll | --lengthlimit
// -fl | --fraglength -> So fragment length in regads to PE and then the current -l, -lf, -ld would be read lengths