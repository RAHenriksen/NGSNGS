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

typedef struct{
  int threads1;
  int threads2;
  size_t reads;
  int Glob_seed;
  const char *OutFormat;
  const char *OutName;
  const char *Reference;
  const char *Seq;
  const char *Adapter1;
  const char *Adapter2;
  const char *QualProfile1;
  const char *QualProfile2;
  const char *Briggs;
  int Length;
  const char *LengthFile;
  const char *Poly;
}argStruct;

float myatof(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse Briggs parameters, provide <nv,Lambda,Delta_s,Delta_d>");
  
  return atof(str);
}

void Sizebreak(char *str){
  fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) \n or parse length distribution file (-lf)");
}

void ErrMsg(int messageno){
  if(messageno == 1){fprintf(stderr,"test");}
  else if(messageno == 2){fprintf(stderr,"test2");}
  else {fprintf(stderr,"test4");}
}


int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 \n\n");
  fprintf(fp,"Usage\n./ngsngs [options] -i <input_reference.fa> -r <Number of reads> -l/-lf <fixed length or length file> -seq <SE/PE> -f <Output format> -o <Output name>\n");
  fprintf(fp,"\nRequired: \n");
  fprintf(fp,"-i | --input: \t\t\t Reference file in .fasta format to sample reads from\n");
  fprintf(fp,"-r | --reads: \t\t\t Number of reads to simulate\n");
  fprintf(fp,"-l | --length: \t\t\t Fixed length of simulated fragments, conflicts with -lf option\n");
  fprintf(fp,"-lf | --lengthfile: \t\t CDF of a length distribution, conflicts with -l option\n");
  fprintf(fp,"-seq | --sequencing: \t\t Simulate single-end or paired-end reads\n");
  fprintf(fp,"\t <SE>\t single-end \n \t <PE>\t paired-end\n");
  fprintf(fp,"-f | --format: \t\t\t File format of the simulated outpur reads\n");
  fprintf(fp,"\t <fa||fasta>\t\t nucletide sequence \n \t <fa.gz||fasta.gz>\t compressed nucletide sequence \n \t <fq||fastq>\t\t nucletide sequence with corresponding quality score \n \t <fq.gz||fastq.gz>\t compressed nucletide sequence with corresponding quality score \n \t <bam>\t\t\t Sequence Alignment Map format\n");
  fprintf(fp,"-o | --output: \t\t\t Prefix of output file name.\n");
  fprintf(fp,"\nOptional: \n");
  fprintf(fp,"-h | --help: \t\t\t Print help page\n");
  fprintf(fp,"-v | --version: \t\t Print help page\n");
  fprintf(fp,"-t1 | --threads1: \t\t Number of threads to use for sampling sequence reads\n");
  fprintf(fp,"-t2 | --threads2: \t\t Number of threads to use write down sampled reads, default = 1\n");
  fprintf(fp,"-s | --seed: \t\t\t Random seed, default value being computer time\n");
  fprintf(fp,"-a1 | --adapter1: \t\t Adapter sequence to add for simulated reads (SE) or first read pair (PE)\n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG \n");
  fprintf(fp,"-a2 | --adapter2: \t\t Adapter sequence to add for second read pair (PE) \n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT \n");
  fprintf(fp,"-p | --poly: \t\t Create Poly(X) tails for reads already containing adapters but below the inferred readcycle length. e.g -p G or -p A \n");
  fprintf(fp,"-q1 | --quality1: \t\t Read Quality profile for simulated reads (SE) or first read pair (PE)\n");
  fprintf(fp,"-q2 | --quality2: \t\t Read Quality profile for for second read pair (PE)\n");
  fprintf(fp,"-b | --briggs: \t\t\t Parameters for the damage patterns using the Briggs model\n");
  fprintf(fp,"\t <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007)\n");
  fprintf(fp,"\t nv: Nick rate pr site \n \t Lambda: Geometric distribution parameter for overhang length\n \t Delta_s: PMD rate in single-strand regions\n \t Delta_s: PMD rate in double-strand regions\n");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->threads1 = 1;
  mypars->threads2 = -1;
  mypars->reads = -1;
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
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  mypars->LengthFile = NULL;
  mypars->Poly = NULL;
  ++argv;
  while(*argv){
    fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->reads = atoi(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->OutName = strdup(*(++argv));
    }
    else if(strcasecmp("-t1",*argv)==0 || strcasecmp("--threads1",*argv)==0){
      mypars->threads1 = atoi(*(++argv));
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      mypars->threads2 = atoi(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv));
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      mypars->Seq = strdup(*(++argv));
      if(strcasecmp("SE",mypars->Seq)!=0 && strcasecmp("PE",mypars->Seq)!=0){HelpPage(stderr);} 
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
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      mypars->OutFormat = strdup(*(++argv));
    }
    else if(strcasecmp("-b",*argv)==0 || strcasecmp("--briggs",*argv)==0){
      mypars->Briggs = strdup(*(++argv)); //double nv, double lambda, double delta_s, double delta -> 0.024,0.36,0.68,0.0097
    }
    else if(strcasecmp("-l",*argv)==0 || strcasecmp("--length",*argv)==0){
      mypars->Length = atoi(*(++argv));
    }
    else if(strcasecmp("-lf",*argv)==0 || strcasecmp("--lengthfile",*argv)==0){
      mypars->LengthFile = strdup(*(++argv));
    }
    else if(strcasecmp("-p",*argv)==0 || strcasecmp("--poly",*argv)==0){
      mypars->Poly = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      exit(0);
    }
    
    // -e1 +2 || --error1 +2 
    // -p || --poly G T
    ++argv;
  }
  return mypars;
}

// ------------------------------ //

int main(int argc,char **argv){
  argStruct *mypars = NULL;
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                        strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else{
    mypars = getpars(argc,argv);
    clock_t t = clock();
    time_t t2 = time(NULL);

    const char *fastafile = mypars->Reference;
    const char* OutputFormat = mypars->OutFormat;
    const char* filename = mypars->OutName; //"chr22_out";
    const char* Seq_Type = mypars->Seq;
    size_t No_reads = mypars->reads;

    if (fastafile == NULL){fprintf(stderr,"Input reference file not recognized, provide -i | --input\n");exit(0);}
    if (OutputFormat == NULL){fprintf(stderr,"Output format not recognized, provide -f | --format : <fa, fa.gz, fq, fq.gz, bam>\n");exit(0);}
    if (filename == NULL){fprintf(stderr,"Output filename not provided, provide -o\n");exit(0);}
    if (Seq_Type == NULL){fprintf(stderr,"Sequence type not provided. provide -seq || --sequence : SE (single-end) or PE (paired-end)\n");exit(0);}
    if (No_reads == -1){fprintf(stderr,"Number of reads to simulate not provided, provide -r \n");exit(0);}
 
    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int Glob_seed = mypars->Glob_seed; 
    int threads1 = mypars->threads1;
    int threads2;
    if (mypars->threads2 == -1){
      threads2 = 1; //threads1
    }
    else{
      threads2 = mypars->threads2;
    }
    
     //1e1;
    
    int Thread_specific_Read = static_cast<int>(No_reads/threads1);

    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",Glob_seed);
    fprintf(stderr,"\t-> Number of threads used for sampling: %d and for writing down %d\n",threads1,threads2);
    fprintf(stderr,"\t-> Number of simulated reads: %zd\n",No_reads);

    const char* Adapt_flag;
    const char* Adapter_1;
    const char* Adapter_2;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      Adapter_2 = mypars->Adapter2;
      //Adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
      //Adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";
      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
    }
    else{
      Adapt_flag = "false";
      if (mypars->Poly != NULL){fprintf(stderr,"Poly tail error: Missing adapter sequence, provide adapter sequence (-a1,-a2) as well\n");exit(0);}
    }

    // QUALITY PROFILES
    const char* QualProfile1; const char* QualProfile2;
    QualProfile1 = mypars->QualProfile1; QualProfile2 = mypars->QualProfile2;
    
    const char* QualStringFlag;
    if (QualProfile1 == NULL){QualStringFlag = "false";}
    else{QualStringFlag = "true";}
    //fprintf(stderr,"qualstring test %s",QualStringFlag);
    if (QualStringFlag == "true"){
      if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0 || strcasecmp("bam",OutputFormat)==0){
        if (Seq_Type == "PE" && QualProfile2 == NULL){
          fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
          exit(0);
        }
      }
    }
    else
    {
      if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){
        fprintf(stderr,"Could not parse the Nucleotide Quality profile(s), for SE provide -q1 for PE provide -q1 and -q2. see helppage (-h). \n");
        exit(0);
      }
    }
        
    int qualstringoffset = 0;
    if(strcasecmp("fq",OutputFormat)==0 || strcasecmp("fq.gz",OutputFormat)==0){qualstringoffset = 33;}

    // LENGTH     
    int FixedSize = mypars->Length;
    const char* Sizefile = mypars->LengthFile;

    if (Sizefile == NULL){
      if (FixedSize == -1){
        fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) or parse length distribution file (-lf) see helppage (-h).\n");
        exit(0);
      }
    }
    if (Sizefile != NULL){
      if (FixedSize != -1){
        fprintf(stderr,"Could not parse both length parameters, provide either fixed length size (-l) or parse length distribution file (-lf) see helppage (-h).\n");
        exit(0);
      }
    }
    
    const char* Briggs_Flag;
    float Param[4];
    if (mypars->Briggs != NULL){
      char* BriggsParam = strdup(mypars->Briggs);
      Param[0] = myatof(strtok(BriggsParam,"\", \t"));
      Param[1] = myatof(strtok(NULL,"\", \t"));
      Param[2] = myatof(strtok(NULL,"\", \t"));
      Param[3] = myatof(strtok(NULL,"\", \t"));
      Briggs_Flag = "True";
    }
    else{Briggs_Flag = "False";}
    
    Create_se_threads(seq_ref,threads1,Glob_seed,Thread_specific_Read,filename,
                      Adapt_flag,Adapter_1,Adapter_2,OutputFormat,Seq_Type,
                      Param,Briggs_Flag,Sizefile,FixedSize,qualstringoffset,
                      QualProfile1,QualProfile2,threads2,QualStringFlag,Polynt);

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

  // OPTIONAL DEALLOCATIONS
  free((char *)mypars->Adapter1);
  free((char *)mypars->Adapter2);
  free((char *)mypars->QualProfile1);
  free((char *)mypars->QualProfile2);
  free((char *)mypars->Briggs);
  delete mypars;

  /*double lambda = 0.36;
  std::default_random_engine generator1(rand()%30000+1); // value between 1 and 30000
  std::geometric_distribution<int> distribution1(lambda); // creates lambda distribution
  std::default_random_engine generator2(rand()%30000+1);
  std::geometric_distribution<int> distribution2(lambda);
  unsigned int gen_1 = (unsigned int) myrandgenmodulo((unsigned int) (1000+2),30000); 
  unsigned int seed = 10;
  for (int i = 0; i < 100000; i++){
    int l = distribution1(generator1);
    int l2 = distribution1(generator2);
    int l3 = Random_geometric_k((int) ((rand_r(&seed)%30000)+1),lambda);
    int l4 = Random_geometric_k((int) ((rand_r(&seed)%30000)+1),lambda);
    fprintf(stderr,"%d \t %d\t %d\t %d\n",l,l2,l3,l4);
  }*/
}


// g++ NGSNGS_func.cpp atomic_fq.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl

//  make HTSSRC=/home/wql443/scratch/htslib/
// ./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -seq PE -f fq -o chr22

//g++ NGSNGS_func.cpp Sampling.cpp ArgParse.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
// ./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 10 -s 1 -t1 1 -lf Size_dist/Size_dist_sampling.txt -b 0.024,0.36,0.68,0.0097 -seq SE -q1 Qual_profiles/AccFreqL150R1.txt -f fq -o chr22
/*
const char* poly_test;
  poly_test = strdup("G");
  char test[1024];
  memset(test,(char) poly_test[0], 150);
  std::cout << test << std::endl;
  char test2[1024] = "TTTTTT";
  strncpy(test, test2, strlen(test2));
  std::cout << test << std::endl;
  std::cout << test2 << std::endl;
  exit(0);
  */