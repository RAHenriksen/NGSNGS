#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"
#include "Sampling.h"
#include "mrand.h"
#include "version.h"
#include <signal.h>
#define LENS 4096
#define MAXBINS 100

int VERBOSE = 1;
int SIG_COND = 1;

typedef struct{
  int threads1; //sampling threads, used internally by pthread_create
  int threads2; //compression threads, used external by bgzf_mt and set_hts_options
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
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 0.5.0 \n\n");
  fprintf(fp,"Usage\n./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf <fixed length or length file> -seq <SE/PE> -f <output format> -o <output name prefix>\n");
  fprintf(fp,"\nExample \n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 2 -s 1 -lf Test_Examples/Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut\n");
  fprintf(fp,"\n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t1 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -q2 Test_Examples/Qual_profiles/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut\n");  
  fprintf(fp,"\n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/DeamSubFile.txt -f fa -o MycoBactFaSEOut\n");    
  fprintf(fp,"\n-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"\nRequired: \n\n");
  fprintf(fp,"-i   | --input: \t\t Reference file in fasta format (.fa,.fasta) to sample reads.\n");
  fprintf(fp,"-r   | --reads: \t\t Number of reads to simulate, conflicts with -c option.\n");
  fprintf(fp,"-c   | --coverage: \t\t Depth of Coverage to simulate, conflics with -r option.\n");
  fprintf(fp,"-l   | --length: \t\t Fixed length of simulated fragments, conflicts with -lf & -ld option.\n");
  fprintf(fp,"-lf  | --lengthfile: \t\t CDF of a length distribution, conflicts with -l & -ld option.\n");
  fprintf(fp,"-ld  | --lengthdist: \t\t Discrete or continuous probability distributions, conflicts with -l & -lf option.\n");
  fprintf(fp,"\teg.\t Uni,40,180 || Norm,80,30 || LogNorm,4,1 || Pois,165 || Exp,0.025 || Gam,20,2\n");
  fprintf(fp,"-seq | --sequencing: \t\t Simulate single-end or paired-end reads.\n");
  fprintf(fp,"\t <SE>\t single-end \n \t <PE>\t paired-end.\n");
  fprintf(fp,"-f   | --format: \t\t File format of the simulated output reads.\n");
  fprintf(fp,"\t <fa||fasta||fa.gz||fasta.gz>\t\t Nucletide sequence w. different compression levels. \n \t <fq||fastq||fq.gz||fastq.gz>\t\t Nucletide sequence with corresponding quality score w. different compression levels. \n \t <sam||bam||cram>\t\t\t Sequence Alignment Map format w. different compression levels.\n");
  fprintf(fp,"-o   | --output: \t\t Prefix of output file name.\n");
  fprintf(fp,"\nOptional: \n");
  fprintf(fp,"\nNucleotide Alterations: \n");
  fprintf(fp,"-bcf: \t\t\t\t Binary Variant Calling Format (.bcf)\n");
  fprintf(fp,"-v:  | --variant: \t\t Specific variants to simulate\n");
  fprintf(fp,"\t eg.\t snp || indel. Default = all\n");
  fprintf(fp,"-b   | --briggs: \t\t Parameters for the damage patterns using the Briggs model.\n");
  fprintf(fp,"\t <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007).\n");
  fprintf(fp,"\t nv: Nick rate pr site. \n \t Lambda: Geometric distribution parameter for overhang length.\n \t Delta_s: PMD rate in single-strand regions.\n \t Delta_d: PMD rate in double-strand regions.\n");
  fprintf(fp,"-mf  | --mismatch: \t\t\t Nucleotide substitution frequency file.\n");
  fprintf(fp,"-ne  | --noerror: \t\t Disabling the nucleotide subsitutions based on nucleotide qualities.\n");
  fprintf(fp,"\nRead Specific: \n");
  fprintf(fp,"-chr | --chromosomes: \t\t Specific chromosomes from input reference file.\n");
  fprintf(fp,"-a1  | --adapter1: \t\t Adapter sequence to add for simulated reads (SE) or first read pair (PE).\n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG \n\n");
  fprintf(fp,"-a2  | --adapter2: \t\t Adapter sequence to add for second read pair (PE). \n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT \n\n");
  fprintf(fp,"-p   | --poly: \t\t\t Create Poly(X) tails for reads, containing adapters with lengths below the inferred readcycle length. \n \t e.g -p G or -p A \n");
  fprintf(fp,"-q1  | --quality1: \t\t Read Quality profile for single-end reads (SE) or first read pair (PE).\n");
  fprintf(fp,"-q2  | --quality2: \t\t Read Quality profile for for second read pair (PE).\n");
  fprintf(fp,"\nOther: \n");
  fprintf(fp,"-t1  | --threads1: \t\t Number of sampling threads, default = 1.\n");
  fprintf(fp,"-t2  | --threads2: \t\t Number of compression threads, default = 0.\n");
  fprintf(fp,"-s   | --seed: \t\t\t Random seed, default = current calendar time (s).\n");
  fprintf(fp,"-rand: \t\t\t\t Pseudo-random number generator, OS specific\n");
  fprintf(fp,"\t e.g. linux || unix -> drand48_r (-rand = 0), not available for MacOS.\n");
  fprintf(fp,"\t APPLE and MacOS (-rand = 1).\n");
  //NoAlign //Indiv
  exit(1);
  return 0;
}

float myatof(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse Briggs parameters, provide <nv,Lambda,Delta_s,Delta_d>");
  
  return atof(str);
}

void Sizebreak(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) \n or parse length distribution file (-lf)");
}

void ErrMsg(double messageno){
  if(messageno == 1.0){fprintf(stderr,"\nInput reference file not recognized, provide -i | --input\n");}
  else if(messageno == 2.0){fprintf(stderr,"\nNumber of reads or depth to simulate not recognized, provide either number of reads (-r) or depth of coverage (-c).\n");}
  else if(messageno == 2.2){fprintf(stderr,"\nUnable to utilize the desired number for depth of coverage (-c), use values above 0.\n");}
  else if(messageno == 2.4){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate (-r), use integers above 0.\n");}
  else if(messageno == 2.6){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate for the provided number of threads, i.e threads > no. reads when no. reads equals 1.\n");}
  else if(messageno == 2.99){fprintf(stderr,"\nCould not parse both the number reads and depth to simulate, provide either number of reads (-r) or depth of coverage (-c).\n");}
  else if(messageno == 3.0){fprintf(stderr,"\nCould not parse the length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 3.2){fprintf(stderr,"\nUnable to simulate reads of the desired fixed length (-l), use integers above 0.\n");}
  else if(messageno == 5.0){fprintf(stderr,"\nCould not parse both length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 6.0){fprintf(stderr,"\nSequence type not provided. provide -seq || --sequence : SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 6.5){fprintf(stderr,"\nSequence type not recognized. provide either SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 7.0){fprintf(stderr,"\nOutput format not recognized, provide -f | --format : <fa, fa.gz, fq, fq.gz, sam, bam, cram>.\n");}
  else if(messageno == 8.0){fprintf(stderr,"\nOutput filename not provided, provide -o.\n");}
  else if(messageno == 9.0){fprintf(stderr,"\nUnable to utilize the provided number of threads, use integers above 0.\n");}
  else if(messageno == 10.0){fprintf(stderr,"\nNucleotide for poly(x) homopolymers not recognized, provide -p : <A,G,C,T,N>.\n");}
  else if(messageno == 11.0){fprintf(stderr,"\nCould not parse the Nucleotide Quality profile(s), for format <fq, fq.gz, sam, bam, cram> provide -q1 for SE and -q1, -q2 for PE.\n");}
  else if(messageno == 12.0){fprintf(stderr,"\nBoth a mismatch file and briggs deamination parameters has been provided. Provide either briggs (-p) og mismatch (-mf).\n");}
  else if(messageno == 13.0){fprintf(stderr,"\nOnly variantion type has been provided (-v). Provide variant calling format (-bcf).\n");}
  else if(messageno == 13.5){fprintf(stderr,"\nProvide variantion type not recognized, provide either snp, indel or all (default)\n");}
  else {fprintf(stderr,"\nError with input parameters, see helppage (-h)");}
  fprintf(stderr,"see helppage (-h)\n");
  exit(0);
}

void WarMsg(double messageno){
  if(messageno == 1.0){fprintf(stderr,"\nWarning: for the output format <fa> the quality profiles (-q1, -q2) remains unused\n");}
  else if(messageno == 2.0){fprintf(stderr,"\nWarning: for the output format <fa, sam, bam, cram> the parameter (-p) is rendered moot without the nucleotide quality profiles (-q1, -q2), since poly(x) tails cannot be added to the reads since the read length cannot be inferred\n");}
  else if(messageno == 3.0){fprintf(stderr,"\nWarning: sequencing errors (-e) are not added to the sequence reads without the nucleotide quality profiles (-q1, -q2) \n");}
  else if(messageno == 4.0){fprintf(stderr,"\nWarning: for the output format <fa> the provided nucleotide qualities (-q1, -q2) will not be used \n");}
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->threads1 = 1;
  mypars->threads2 = -1;
  mypars->nreads = 0;
  mypars->coverage = -1.0;
  mypars->Length = -1;
  mypars->Glob_seed = (int) time(NULL);
  mypars->OutFormat = NULL; //"fa";
  mypars->OutName = NULL; //"output";
  mypars->Seq = NULL; // "SE";
  mypars->Reference = NULL;
  mypars->Adapter1 = NULL;
  mypars->Adapter2 = NULL;
  mypars->QualProfile1 = NULL;
  mypars->QualProfile2 = NULL;
  mypars->SubProfile = NULL;
  mypars->ErrorFlag = NULL;
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  mypars->LengthFile = NULL;
  mypars->LengthDist = NULL;
  mypars->Chromosomes = NULL;
  mypars->Poly = NULL;
  mypars->rand_val = -1;
  mypars->Variant = NULL;
  mypars->Variant_type = NULL;
  mypars->HeaderIndiv=NULL;
  mypars->NoAlign=NULL;
  mypars->CommandRun = (char*) calloc(1024,1);
  char *Command = mypars->CommandRun;
  const char *first = "./ngsngs ";
  strcpy(Command,first);
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
      strcat(Command,mypars->Variant); strcat(Command," ");
      if(mypars->Variant == NULL){ErrMsg(13.0);}
      else if (mypars->Variant_type && strcasecmp("snp",mypars->Variant_type)!=0 && 
      strcasecmp("indel",mypars->Variant_type)!=0 &&
      strcasecmp("all",mypars->Variant_type)!=0){ErrMsg(13.5); exit(0);}      
    }
    else if(strcasecmp("-t1",*argv)==0 || strcasecmp("--threads1",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->threads1 = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->threads1 < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->threads2 = atoi(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->threads2 < 0){ErrMsg(9.0);}
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      const char* readstr = strdup(*(++argv));
      sscanf(readstr, "%lu",&mypars->nreads);
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->nreads <= 0){ErrMsg(2.4);}
    }
    else if(strcasecmp("-c",*argv)==0 || strcasecmp("--cov",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->coverage = atof(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
      if (mypars->coverage <= 0.0){ErrMsg(2.2);}
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
      mypars->ErrorFlag = "F"; // strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-na",*argv)==0 || strcasecmp("--noalign",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->NoAlign = "T"; // strdup(*(++argv));
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
      if (mypars->Length <= 0.0){ErrMsg(3.2);}
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
    else if(strcasecmp("-indiv",*argv)==0){
      strcat(Command,*argv); strcat(Command," ");
      mypars->HeaderIndiv = strdup(*(++argv));
      strcat(Command,*argv); strcat(Command," ");
    }
    else if(strcasecmp("-rand",*argv)==0){
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
  return mypars;
}


int really_kill =3;

void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
    fprintf(stderr,"\n\t-> If you really want ngsngs to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;
  //  pthread_mutex_unlock(&mUpPile_mutex);
}

void catchkill(){
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  
}


// ------------------------------ //

int main(int argc,char **argv){
  //fprintf(stderr,"PRINT TYPE %d\n",MacroRandType);
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    catchkill();
    mypars = getpars(argc,argv);
    if(mypars==NULL)
      return 1;
    
    char* Command = mypars->CommandRun;
    fprintf(stderr,"\n\t-> ngsngs version: %s (htslib: %s) build(%s %s)\n",NGSNGS_VERSION,hts_version(),__DATE__,__TIME__); 
    fprintf(stderr,"\t-> Mycommmand: %s\n strlen: %lu\n",Command,strlen(mypars->CommandRun));

    //fprintf(stderr,"\t-> Command 2 : %s and version %s \n",CommandArray,version);
    clock_t t = clock();
    time_t t2 = time(NULL);
    int Glob_seed = mypars->Glob_seed; 

    const char *fastafile = mypars->Reference;
    const char* OutputFormat = mypars->OutFormat;
    const char* filename = mypars->OutName; //"chr22_out";
    const char* Seq_Type = mypars->Seq;
    //    size No_reads = mypars->nreads;
    double readcov = mypars->coverage;
    int MacroRandType = mypars->rand_val; //extern int

    if (MacroRandType == -1){
      #if defined(__linux__) || defined(__unix__) // all unices not caught above
      // Unix
        MacroRandType = 0;
      #elif defined(__APPLE__) || defined(__MACH__)
        MacroRandType = 1;
      #else
      #   error "Unknown compiler"
      #endif
    }
    //fprintf(stderr,"RANDOM VALUE %d \n",MacroRandType);
    
    
    if (fastafile == NULL){ErrMsg(1.0);}
    if (Seq_Type == NULL){ErrMsg(6.0);}
    if (OutputFormat == NULL){ErrMsg(7.0);}
    if (filename == NULL){ErrMsg(8.0);}
    
    if (mypars->nreads == 0 && readcov== -1.0 )
        ErrMsg(2.0);

    //fprintf(stderr,"\t-> Command: %s \n",Command);
    int FixedSize = mypars->Length;
    const char* Sizefile = mypars->LengthFile;
    const char* SizeDist = mypars->LengthDist;
    double meanlength = 0;// DRAGEON REMEMBEWR TO CHECK meanlength is interpreted as float
    int SizeDistType=-1;int val1; int val2;

    if (FixedSize != -1){
      if (FixedSize == -1){ErrMsg(3.0);}
      else{
        meanlength = FixedSize;
      } 
    }
    if (Sizefile != NULL){
      //fprintf(stderr,"SIZE FILE ARG\n");
      double sum,n;
      sum=n=0;

      char buf[LENS];
      gzFile gz = Z_NULL;
      gz = gzopen(Sizefile,"r");
      assert(gz!=Z_NULL);
      while(gzgets(gz,buf,LENS)){
        double Length_tmp = atof(strtok(buf,"\n\t "));
        double Frequency_tmp = atof(strtok(NULL,"\n\t "));
        sum += Length_tmp*Frequency_tmp;
        n = n+1;
      }
      gzclose(gz);
      
      meanlength = sum/n; //mindfuck
      if (FixedSize != -1){ErrMsg(5.0);}
    }

    
    if (SizeDist != NULL){
      //fprintf(stderr,"LENGTH DISTRIBUTION and %s",SizeDist);
      std::default_random_engine generator(Glob_seed);
      char* Dist;
      
      char* DistParam = strdup(SizeDist);
      Dist = strtok(DistParam,",");
      val1 = atoi(strtok (NULL, ","));
      //fprintf(stderr,"strtok %d\n",val1);
      char* tmp = strtok(NULL, ",");
      if(tmp == NULL){val2 = 0;}
      else{val2 = atoi(tmp);}
      
      if (strcasecmp(Dist,"Uni")==0){SizeDistType=1;std::uniform_int_distribution<int> distribution(val1,val2);meanlength=(0.5*(val1+val2));}
      if (strcasecmp(Dist,"Norm")==0){SizeDistType=2;std::normal_distribution<double> distribution(val1,val2);meanlength= val1;}
      if (strcasecmp(Dist,"LogNorm")==0){SizeDistType=3;std::lognormal_distribution<double> distribution(val1,val2);meanlength= exp((val1+((val2*val2)/2)));}
      if (strcasecmp(Dist,"Pois")==0){SizeDistType=4;std::poisson_distribution<int> distribution(val1);meanlength= val1;}
      if (strcasecmp(Dist,"Exp")==0){SizeDistType=5;std::exponential_distribution<double> distribution(val1);meanlength= 1/val1;}
      if (strcasecmp(Dist,"Gam")==0){SizeDistType=6;std::gamma_distribution<double> distribution(val1,val2);meanlength= (val1/val2);}
      if (FixedSize != -1){ErrMsg(5.0);}
      free((char *)Dist);
    }

    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int threads1 = mypars->threads1;
    int threads2;
    if (mypars->threads2 == -1){threads2 = 1;}
    else{threads2 = mypars->threads2;}


    //first capture the awkward cases where no reads or cov has been defined or both defined
    if(mypars->nreads == 0 &&readcov ==-1.0){
      fprintf(stderr,"must suply nreads or cov");exit(0);
    }
    if(mypars->nreads > 0 &&readcov !=-1.0){
      fprintf(stderr,"must not suply nreads and cov");exit(0);
    }
    //now compute the number of reads required acroos all threads

    if (readcov != -1){
      size_t genome_size = 0;

      for (int i = 0; i < chr_total; i++){
        const char *chr_name = faidx_iseq(seq_ref,i);
        int chr_len = faidx_seq_len(seq_ref,chr_name);
        genome_size += chr_len;
      }
      mypars->nreads =  (readcov*genome_size)/meanlength;
      fprintf(stderr,"\t-> Number of simulated reads: %zu or coverage: %f\n",mypars->nreads,mypars->coverage);
    }

    size_t nreads_per_thread = mypars->nreads/threads1;
    
    //fprintf(stderr,"\t-> Command: %s \n",Command);
    //fprintf(stderr,"\t-> Command 2 : %s \n",CommandArray);
    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",Glob_seed);
    fprintf(stderr,"\t-> Number of sampling threads used (-t1): %d and number of compression threads (-t2): %d\n",threads1,threads2);
    fprintf(stderr,"\t-> Number of simulated reads: %zu or coverage: %f\n",mypars->nreads,mypars->coverage);

    const char* Adapt_flag;
    const char* Adapter_1;
    const char* Adapter_2;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      //fprintf(stderr,"\t-> ARGPARSE ADAPTER + POLY\n");
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      Adapter_2 = mypars->Adapter2;

      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
      else{Polynt = "F";}
      
    }
    else{
      //fprintf(stderr,"\t-> ARGPARSE ADAPT FLAG+ POLY\n");
      Adapt_flag = "false";
      if (mypars->Poly != NULL){fprintf(stderr,"Poly tail error: Missing adapter sequence, provide adapter sequence (-a1,-a2) as well\n");exit(0);}
      else{Polynt = "F";}
    }
    // QUALITY PROFILES
    const char* QualProfile1; const char* QualProfile2;
    QualProfile1 = mypars->QualProfile1; QualProfile2 = mypars->QualProfile2;

    const char* QualStringFlag;
    if (QualProfile1 == NULL){QualStringFlag = "false";}
    else{QualStringFlag = "true";}
    //fprintf(stderr,"qualstring test %s",QualStringFlag);
    if (strcasecmp("true",QualStringFlag)==0){
      if(OutputFormat && (strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0 || strcasecmp("sam",OutputFormat)==0 || strcasecmp("bam",OutputFormat)==0 || strcasecmp("cram",OutputFormat)==0)){
        if (strcasecmp("PE",Seq_Type)==0 && QualProfile2 == NULL){
          ErrMsg(11.0);
          //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
          exit(0);
        }
      }
    }
    else
    {
      if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){
        ErrMsg(11.0);
        //fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
    }
    //fprintf(stderr,"\t-> ADAPTER FLAG IS:%s\n",Adapt_flag);
    //fprintf(stderr,"\t-> QUAL STRING FLAG IS:%s\n",QualStringFlag);
    const char* ErrorFlag;
    if (mypars->ErrorFlag != NULL){
      ErrorFlag = mypars->ErrorFlag;
    }
    else{ErrorFlag = "T";}
    
    const char* NoAlign;
    if (mypars->NoAlign != NULL){
      NoAlign = mypars->NoAlign;
    }
    else{NoAlign = "F";}

    if (strcasecmp("false",QualStringFlag)==0){
      if(strcasecmp("true",Adapt_flag)==0 && mypars->Poly != NULL){WarMsg(2.0);}
      if(mypars->ErrorFlag != NULL){WarMsg(3.0);}
    }
    if (strcasecmp("true",QualStringFlag)==0){
      if(strcasecmp("fa",OutputFormat)==0 || strcasecmp("fa.gz",OutputFormat)==0){WarMsg(4.0);}
    }

    int qualstringoffset = 0;
    if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){qualstringoffset = 33;}
    
    const char* Briggs_Flag;
    float Param[4];
    if (mypars->Briggs != NULL){
      char* BriggsParam = strdup(mypars->Briggs);
      Param[0] = myatof(strtok(BriggsParam,"\", \t"));
      Param[1] = myatof(strtok(NULL,"\", \t"));
      Param[2] = myatof(strtok(NULL,"\", \t"));
      Param[3] = myatof(strtok(NULL,"\", \t"));
      Briggs_Flag = "True";
      free(BriggsParam); // Again using strdup
    }
    else{Briggs_Flag = "False";}
    
    const char* SubProfile; const char* SubFlag;
    
    SubProfile = mypars->SubProfile;
    if (SubProfile == NULL){SubFlag = "false";}
    else{SubFlag = "true";}
    //fprintf(stderr,"SUB FLAG IS %s\n",SubFlag);
    if(SubProfile != NULL && mypars->Briggs != NULL){
      ErrMsg(12.0);
      exit(0);
    }

    const char* Specific_Chr[1024] = {};
    if (mypars->Chromosomes != NULL){
      fprintf(stderr,"PARTIAL chromosomes %s\n",mypars->Chromosomes);
      int chr_idx_partial = 0;
      Specific_Chr[chr_idx_partial++] = strtok(strdup(mypars->Chromosomes),"\", \t");
      char *chrtok = NULL;
      while(((chrtok=strtok(NULL,"\", \t")))){
	      Specific_Chr[chr_idx_partial++] = strdup(chrtok);
	      assert(chr_idx_partial<MAXBINS);
      }
      fprintf(stderr,"AFTER WHILE and chr_idx_partial %d\n",chr_idx_partial);
      Specific_Chr[chr_idx_partial++] = "\0";
    }
    
    char* Variant_flag =NULL;
    const char* VCFformat = mypars->Variant;
    const char* VarType =mypars->Variant_type;
    const char* HeaderIndiv = mypars->HeaderIndiv;

    if (VCFformat != NULL){
      Variant_flag = strdup("bcf");
      if(VarType == NULL){
        VarType = strdup("snp");
      }
      else if(VarType != NULL){
        VarType = mypars->Variant_type;
      }
      fprintf(stderr,"VARIANT TYPE %s\n",VarType);
    }

    //if(Specific_Chr[0]=='\0'){fprintf(stderr,"HURRA");}
    int DeamLength = 0;
    //const char* HeaderIndiv = "HG00096";
    Create_se_threads(seq_ref,threads1,Glob_seed,nreads_per_thread,filename,
                      Adapt_flag,Adapter_1,Adapter_2,OutputFormat,Seq_Type,
                      Param,Briggs_Flag,Sizefile,FixedSize,SizeDistType,val1,val2,
                      qualstringoffset,QualProfile1,QualProfile2,threads2,QualStringFlag,Polynt,
                      ErrorFlag,Specific_Chr,fastafile,SubFlag,SubProfile,DeamLength,MacroRandType,
                      VCFformat,Variant_flag,VarType,Command,NGSNGS_VERSION,HeaderIndiv,NoAlign);
    fai_destroy(seq_ref); //ERROR SUMMARY: 8 errors from 8 contexts (suppressed: 0 from 0) definitely lost: 120 bytes in 5 blocks
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));
  }

  // MEMORY DEALLOCATION OF STRDUP FROM INPUT PARAMETERS
  // REQUIRED DEALLOCATIONS
  free((char *)mypars->Reference);
  free((char *)mypars->Seq);
  
  free((char *)mypars->OutFormat);
  free((char *)mypars->OutName);
  free((char *)mypars->LengthFile);
  free((char *)mypars->LengthDist);
  
  // OPTIONAL DEALLOCATIONS
  free((char *)mypars->Adapter1);
  free((char *)mypars->Adapter2);
  free((char *)mypars->QualProfile1);
  free((char *)mypars->QualProfile2);
  free((char *)mypars->SubProfile);
  free((char *)mypars->Briggs);
  free((char *)mypars->Poly);
  free((char *)mypars->HeaderIndiv);
  delete mypars;
}
