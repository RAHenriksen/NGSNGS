#include <htslib/kstring.h>
#include "NGSNGS_cli.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include "version.h"
#include <cassert>

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->SamplThreads = 1;
  mypars->CompressThreads = 1;
  int compress_t_y_n = 0; // set to no pr default

  // Simulation mode
  mypars->simmode = 0;

  // generating strings for which the simulated reads will be contained
  mypars->nreads = 0;
  mypars->coverage = 0.0;
  mypars->KstrBuf=30000000;

  // The output format, output files, and structural elements for SAM outputs
  mypars->OutFormat = unknownT;
  mypars->OutName = NULL;
  mypars->DumpFile = NULL;
  mypars->IndelDumpFile = NULL;
  mypars->HeaderIndiv=-1;
  mypars->NameIndiv = NULL;
  mypars->Align=1;

  // Thread generation and sampling specific information
  mypars->Chromosomes = NULL;
  mypars->Reference = NULL;
  mypars->seq_type = unknownTT;
  mypars->Glob_seed = (int) time(NULL);
  mypars->Glob_seed_binary = 0;
  mypars->rng_type = -1;
  mypars->BedFile = NULL;
  mypars->flankingregion = 30;
  mypars->MaskBed = 0;
  mypars->CaptureVCF = 0;
  mypars->linkage = 0;

  // Sequence alteration models
  // 1) nucleotide quality score and sequencing errors,  
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  mypars->FixedQual = 0;
  mypars->DoSeqErr = 1;

  // 2) briggs model
  mypars->Briggs = NULL;
  mypars->BriggsBiotin = NULL;
  mypars->Duplicates = 1;
  
  // 3) misincorporation matrix
  mypars->SubProfile = NULL;
  mypars->MisMatchMatrix_bdam = NULL;
  mypars->M3outname = NULL;
  // 4) Bcf file and variation incorporation
  mypars->vcffile = NULL;

  // 5) random variations to the reference genome
  mypars->mutationrate = 0.0;
  mypars->generations = 1;
  mypars->referencevariations = 0;

  // Fragment lengths
  mypars->CycleLength = 0; 
  mypars->LowerLimit=30;
  mypars->Length = 0;
  mypars->LengthFile = NULL;
  mypars->LengthDist = NULL;

  // Additional information for sequence reads
  mypars->Adapter1 = NULL;
  mypars->Adapter2 = NULL;
  mypars->Poly = NULL;

  mypars->Indel = NULL;

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
      if(strcasecmp("cram",tok)==0)
	      mypars->OutFormat = cramT;
      if(mypars->OutFormat==unknownT){
	      fprintf(stderr,"\nNext Generation Simulator for Next Generator Sequencing Data\nWarning:\n");
        ErrMsg(7.0);
      }
    }
    else if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-vcf",*argv)==0 || strcasecmp("-bcf",*argv)==0){
      mypars->vcffile = strdup(*(++argv));
    }
    else if(strcasecmp("-t",*argv)==0 || strcasecmp("--threads",*argv)==0){
      mypars->SamplThreads = atoi(*(++argv));
      if (mypars->SamplThreads < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      compress_t_y_n = 1; //yes the compression are using the threads
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
    else if(strcasecmp("-c",*argv)==0 || strcasecmp("--coverage",*argv)==0){
      mypars->coverage = atof(*(++argv));
      if (mypars->coverage <= 0.0){ErrMsg(2.2);}
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      /*if(*(++argv) == NULL){
 	      fprintf(stderr,"\nNext Generation Simulator for Next Generator Sequencing Data\nWarning:\n");
        ErrMsg(8.1);
        exit(1);
      }*/
      mypars->OutName = strdup(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv))*1000;
      mypars->Glob_seed_binary = 1;
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      char * tok = *(++argv);
      
      if(strcasecmp("SE",tok)==0 || strcasecmp("se",tok)==0 || strcasecmp("single",tok)==0 || strcasecmp("single-end",tok)==0)
	      mypars->seq_type = SE;
      else if(strcasecmp("PE",tok)==0 || strcasecmp("pe",tok)==0 || strcasecmp("paired",tok)==0 || strcasecmp("paired-end",tok)==0)
	      mypars->seq_type = PE;
      else if(mypars->seq_type==unknownTT)
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
    else if(strcasecmp("-qs",*argv)==0 || strcasecmp("--qualityscore",*argv)==0){
      mypars->FixedQual = atoi(*(++argv));
    }
    else if(strcasecmp("-mf",*argv)==0 || strcasecmp("--mismatch",*argv)==0){
      mypars->SubProfile = strdup(*(++argv));
    }
    else if(strcasecmp("-m3",*argv)==0 || strcasecmp("--mismatchmatrix",*argv)==0){
      mypars->MisMatchMatrix_bdam = strdup(*(++argv));
    }
    else if(strcasecmp("-Dumpm3",*argv)==0){
      mypars->M3outname = strdup(*(++argv));
    }
    else if(strcasecmp("-ne",*argv)==0 || strcasecmp("--noerror",*argv)==0){
      mypars->DoSeqErr = 0;//"F";
    }
    else if(strcasecmp("-na",*argv)==0 || strcasecmp("--noalign",*argv)==0){
      mypars->Align = 0;
    }
    else if(strcasecmp("-m",*argv)==0 || strcasecmp("--model",*argv)==0){
      ++argv;
      char *tok = *argv;
      char* ModelString = strdup(tok);
      char* BriggsModel;
      BriggsModel = strtok(ModelString,",");
      char* ModelParam =  strdup(strtok (NULL, ""));
      if(strcasecmp("b",BriggsModel)==0 || strcasecmp("briggs",BriggsModel)==0){
	      mypars->Briggs = ModelParam;
      }
      if(strcasecmp("b7",BriggsModel)==0 || strcasecmp("briggs07",BriggsModel)==0)
	      mypars->BriggsBiotin = ModelParam;
      free(ModelString);
    }
    else if(strcasecmp("-dup",*argv)==0 || strcasecmp("--duplicates",*argv)==0){
      mypars->Duplicates = atoi(*(++argv));
    }
    else if(strcasecmp("-cl",*argv)==0 || strcasecmp("--cycle",*argv)==0){
      mypars->CycleLength = atoi(*(++argv));
      if (mypars->CycleLength < 0.0){ErrMsg(3.2);}
    }
    else if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      mypars->Length = atoi(*(++argv));
      if (mypars->Length < 0.0){ErrMsg(3.1);}
    }
    else if(strcasecmp("-lf",*argv)==0 || strcasecmp("--lengthfile",*argv)==0){
      mypars->LengthFile = strdup(*(++argv));
    }
    else if(strcasecmp("-ld",*argv)==0 || strcasecmp("--lengthdist",*argv)==0){
      mypars->LengthDist = strdup(*(++argv));
    }
    else if(strcasecmp("-ll",*argv)==0 || strcasecmp("--lowerlimit",*argv)==0){
      mypars->LowerLimit = atoi(*(++argv));
      if (mypars->LowerLimit < 30.0){ErrMsg(3.3);}
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
      mypars->HeaderIndiv = atoi(*(++argv));
    }
    else if(strcasecmp("-name",*argv)==0 || strcasecmp("--headername",*argv)==0){
      mypars->NameIndiv = strdup(*(++argv));;
    }
    else if(strcasecmp("-rng",*argv)==0 || strcasecmp("--rand",*argv)==0){
      mypars->rng_type = atoi(*(++argv));
    }
    else if(strcasecmp("-indel",*argv)==0){
      mypars->Indel = strdup(*(++argv));
    }
    else if(strcasecmp("-DumpVCF",*argv)==0){
      mypars->DumpFile = strdup(*(++argv));
    }
    else if(strcasecmp("-DumpIndel",*argv)==0){
      mypars->IndelDumpFile = strdup(*(++argv));
    }
    else if(strcasecmp("-mr",*argv)==0 || strcasecmp("--mutationrate",*argv)==0){
      mypars->mutationrate = atof(*(++argv));
    }
    else if(strcasecmp("-g",*argv)==0 || strcasecmp("--generations",*argv)==0){
      mypars->generations = atoi(*(++argv));
    }
    else if(strcasecmp("-v",*argv)==0 || strcasecmp("--variations",*argv)==0){
      mypars->referencevariations = atol(*(++argv));
    }
    else if(strcasecmp("-circ",*argv)==0 || strcasecmp("--circular",*argv)==0){
      mypars->simmode = 1;
    }
    else if(strcasecmp("-fl",*argv)==0 || strcasecmp("--flanking",*argv)==0){
      mypars->flankingregion = atol(*(++argv));
    }
    else if(strcasecmp("-incl",*argv)==0 || strcasecmp("--include",*argv)==0){
      mypars->BedFile = strdup(*(++argv));
    }
    else if(strcasecmp("-excl",*argv)==0 || strcasecmp("--exclude",*argv)==0){
      mypars->BedFile = strdup(*(++argv));
      mypars->MaskBed = 1;
    }
    else if(strcasecmp("-cap",*argv)==0 || strcasecmp("--capture",*argv)==0){
      mypars->CaptureVCF = 1;
    }
    else if(strcasecmp("-linkage",*argv)==0 || strcasecmp("--linkagedisequilibrium",*argv)==0){
      mypars->linkage = 1;
    }
    else{
      fprintf(stderr,"Unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      return NULL;
    }
     
    ++argv;
  }

  // adjust the compression threads following the input parameters depending on the sampling threads and output file format
  if (mypars->OutFormat == fagzT || mypars->OutFormat == fqgzT || mypars->OutFormat == bamT || mypars->OutFormat == cramT){
    if (mypars->SamplThreads <= 12 && compress_t_y_n == 0){
        mypars->CompressThreads = mypars->SamplThreads;
    }
    else if (mypars->SamplThreads > 12 && compress_t_y_n == 0){
      //putting an upper limit on the number of compression threads
      mypars->CompressThreads = 12;
    }
  }

  if(mypars->linkage == 1 && mypars->CaptureVCF == 1){
    fprintf(stderr,"Conflicting simulations for genome-of-interest with variants, please provide either -linkage or -capture");
    exit(1);  
  }

  // Ensure the required argument are parsed after all arguments have been provided

  const char* NGSNGS_msg = "Next Generation Simulator for Next Generator Sequencing Data";
  // Input
  if(mypars->Reference == NULL){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(1.0);
    exit(1);
  }
  // read numbers
  if(mypars->nreads == 0 && mypars->coverage == 0.0){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(2.5);
  }
  // read lengths
  if(mypars->Length == 0 && mypars->LengthFile == NULL && mypars->LengthDist == NULL){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(3.0);
  }
  // Format
  if(mypars->OutFormat == unknownT){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(7.0);
   exit(1);
  }
  if(mypars->seq_type == unknownTT){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(6.0);
    exit(1);
  }
  // quality profiles
  if(mypars->OutFormat==fqT|| mypars->OutFormat== fqgzT ||mypars->OutFormat==samT ||mypars->OutFormat==bamT|| mypars->OutFormat== cramT){
    if (mypars->QualProfile1 == NULL && mypars->FixedQual == 0){
      fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
      ErrMsg(11.0);
      exit(1);
    }
  }
  // Output
  if(mypars->OutName == NULL){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(8.0);
    exit(1);
  }
  //reference variations
  if(mypars->mutationrate>0.0 && mypars->referencevariations>0){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(15.0);
    //if(mypars->mutationrate<0.0 || mypars->referencevariations<0){ErrMsg(15.1);}
  }

  //vcf
  if(mypars->HeaderIndiv>=0 && mypars->NameIndiv!=NULL){
    fprintf(stderr,"\n%s\nWarning:\n",NGSNGS_msg);
    ErrMsg(16.0);
  }
  
  return mypars;
}

void argStruct_destroy(argStruct *mypars){
  free(mypars->Reference);
  free(mypars->OutName);
  free(mypars->DumpFile);
  free(mypars->IndelDumpFile);
  free(mypars->LengthFile);
  free(mypars->LengthDist);
  free(mypars->Indel);
  
  free(mypars->CommandRun);
  if(mypars->Chromosomes)
    free(mypars->Chromosomes);
  if(mypars->BedFile)
    free(mypars->BedFile);
  free(mypars->vcffile);
  free(mypars->Adapter1);
  free(mypars->Adapter2);
  free(mypars->QualProfile1);
  free(mypars->QualProfile2);
  free(mypars->SubProfile);
  free(mypars->MisMatchMatrix_bdam);
  free(mypars->M3outname);
  free(mypars->Briggs);
  free(mypars->BriggsBiotin);
  free(mypars->Poly);
  delete mypars;
}