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
#define LENS 4096
#define MAXBINS 100

//#ifdef __APPLE__

#if defined(__APPLE__) && defined(__MACH__) 
#include "NGSNGS_Random.h"
#endif /* __APPLE__ */

typedef struct{
  int threads1;
  int threads2;
  size_t reads;
  double coverage;
  int Glob_seed;
  const char *OutFormat;
  const char *OutName;
  const char *Reference;
  const char *Seq;
  const char *Adapter1;
  const char *Adapter2;
  const char *QualProfile1;
  const char *QualProfile2;
  const char *ErrorFlag;
  const char *Briggs;
  int Length;
  const char *LengthFile;
  const char *Poly;
  const char *Chromosomes;
}argStruct;


int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 \n\n");
  fprintf(fp,"Usage\n./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf <fixed length or length file> -seq <SE/PE> -f <output format> -o <output name prefix>\n");
  fprintf(fp,"\n Options: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-i   | --input: \t\t Reference file in fasta format (.fa,.fasta) to sample reads.\n");
  fprintf(fp,"-chr | --chromosomes: \t\t Specific chromosomes from input reference file.\n");
  fprintf(fp,"-r   | --reads: \t\t Number of reads to simulate, conflicts with -c option.\n");
  fprintf(fp,"-c   | --coverage: \t\t Depth of Coverage to simulate, conflics with -r option.\n");
  fprintf(fp,"-l   | --length: \t\t Fixed length of simulated fragments, conflicts with -lf option.\n");
  fprintf(fp,"-lf  | --lengthfile: \t\t CDF of a length distribution, conflicts with -l option.\n");
  fprintf(fp,"-seq | --sequencing: \t\t Simulate single-end or paired-end reads.\n");
  fprintf(fp,"\t <SE>\t single-end \n \t <PE>\t paired-end.\n");
  fprintf(fp,"-f   | --format: \t\t File format of the simulated output reads.\n");
  fprintf(fp,"\t <fa||fasta>\t\t Nucletide sequence. \n \t <fa.gz||fasta.gz>\t Compressed nucletide sequence. \n \t <fq||fastq>\t\t Nucletide sequence with corresponding quality score. \n \t <fq.gz||fastq.gz>\t Compressed nucletide sequence with corresponding quality score. \n \t <bam>\t\t\t Sequence Alignment Map format.\n");
  fprintf(fp,"-o   | --output: \t\t Prefix of output file name.\n");
  fprintf(fp,"-t1  | --threads1: \t\t Number of threads to use for sampling sequence reads.\n");
  fprintf(fp,"-t2  | --threads2: \t\t Number of threads to use write down sampled reads, default = 1.\n");
  fprintf(fp,"-s   | --seed: \t\t\t Random seed, default = current calendar time (s).\n");
  fprintf(fp,"-a1  | --adapter1: \t\t Adapter sequence to add for simulated reads (SE) or first read pair (PE).\n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG \n\n");
  fprintf(fp,"-a2  | --adapter2: \t\t Adapter sequence to add for second read pair (PE). \n");
  fprintf(fp,"\t e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT \n\n");
  fprintf(fp,"-p   | --poly: \t\t\t Create Poly(X) tails for reads, containing adapters with lengths below the inferred readcycle length. \n \t e.g -p G or -p A \n");
  fprintf(fp,"-q1  | --quality1: \t\t Read Quality profile for single-end reads (SE) or first read pair (PE).\n");
  fprintf(fp,"-q2  | --quality2: \t\t Read Quality profile for for second read pair (PE).\n");
  fprintf(fp,"-b   | --briggs: \t\t Parameters for the damage patterns using the Briggs model.\n");
  fprintf(fp,"\t <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007).\n");
  fprintf(fp,"\t nv: Nick rate pr site. \n \t Lambda: Geometric distribution parameter for overhang length.\n \t Delta_s: PMD rate in single-strand regions.\n \t Delta_s: PMD rate in double-strand regions.\n");
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
  else if(messageno == 2.2){fprintf(stderr,"\nUnable to utilize the desired number for depth of coverage, use values above 0.\n");}
  else if(messageno == 2.4){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate, use integers above 0.\n");}
  else if(messageno == 2.6){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate for the provided number of threads, i.e threads > no. reads when no. reads equals 1.\n");}
  else if(messageno == 3.0){fprintf(stderr,"\nCould not parse both the number reads and depth to simulate, provide either number of reads (-r) or depth of coverage (-c).\n");}
  else if(messageno == 4.0){fprintf(stderr,"\nCould not parse the length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 5.0){fprintf(stderr,"\nCould not parse both length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 6.0){fprintf(stderr,"\nSequence type not provided. provide -seq || --sequence : SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 6.5){fprintf(stderr,"\nSequence type not recognized. provide either SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 7.0){fprintf(stderr,"\nOutput format not recognized, provide -f | --format : <fa, fa.gz, fq, fq.gz, sam, bam>.\n");}
  else if(messageno == 8.0){fprintf(stderr,"\nOutput filename not provided, provide -o.\n");}
  else if(messageno == 9.0){fprintf(stderr,"\nUnable to utilize the provided number of threads, use integers above 0.\n");}
  else if(messageno == 10.0){fprintf(stderr,"\nNucleotide for poly(x) homopolymers not recognized, provide -p : <A,G,C,T,N>.\n");}
  else {fprintf(stderr,"\nError with input parameters, see helppage (-h)");}
  fprintf(stderr,"see helppage (-h)\n");
  exit(0);
}



argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->threads1 = 1;
  mypars->threads2 = -1;
  mypars->reads = -1;
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
  mypars->ErrorFlag = NULL;
  mypars->Briggs = NULL; //"0.024,0.36,0.68,0.0097";
  mypars->LengthFile = NULL;
  mypars->Chromosomes = NULL;
  mypars->Poly = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->Reference = strdup(*(++argv));
    }
    else if(strcasecmp("-t1",*argv)==0 || strcasecmp("--threads1",*argv)==0){
      mypars->threads1 = atoi(*(++argv));
      if (mypars->threads1 < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-t2",*argv)==0 || strcasecmp("--threads2",*argv)==0){
      mypars->threads2 = atoi(*(++argv));
      if (mypars->threads2 < 1){ErrMsg(9.0);}
    }
    else if(strcasecmp("-r",*argv)==0 || strcasecmp("--reads",*argv)==0){
      mypars->reads = atoi(*(++argv));
      if (mypars->reads <= 0){ErrMsg(2.4);}
    }
    else if(strcasecmp("-c",*argv)==0 || strcasecmp("--cov",*argv)==0){
      mypars->coverage = atof(*(++argv));
      if (mypars->coverage <= 0.0){ErrMsg(2.2);}
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->OutName = strdup(*(++argv));
    }
    else if(strcasecmp("-s",*argv)==0 || strcasecmp("--seed",*argv)==0){
      mypars->Glob_seed = atoi(*(++argv));
    }
    else if(strcasecmp("-seq",*argv)==0 || strcasecmp("--sequencing",*argv)==0){
      mypars->Seq = strdup(*(++argv));
      if(strcasecmp("SE",mypars->Seq)!=0 && strcasecmp("PE",mypars->Seq)!=0){ErrMsg(6.5);} 
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
    else if(strcasecmp("-e",*argv)==0 || strcasecmp("--error",*argv)==0){
      mypars->ErrorFlag = "T"; // strdup(*(++argv));
    }
    else if(strcasecmp("-f",*argv)==0 || strcasecmp("--format",*argv)==0){
      mypars->OutFormat = strdup(*(++argv));
      if(strcasecmp("fa",mypars->OutFormat)!=0 && 
      strcasecmp("fa.gz",mypars->OutFormat)!=0 &&
      strcasecmp("fq",mypars->OutFormat)!=0 &&
      strcasecmp("fq.gz",mypars->OutFormat)!=0 &&
      strcasecmp("sam",mypars->OutFormat)!=0 &&
      strcasecmp("bam",mypars->OutFormat)!=0){ErrMsg(7.0);}      
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
    else{
      fprintf(stderr,"unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      exit(0);
    }
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
    double readcov = mypars->coverage;

    if (fastafile == NULL){ErrMsg(1.0);}
    if (Seq_Type == NULL){ErrMsg(6.0);}
    if (OutputFormat == NULL){ErrMsg(7.0);}
    if (filename == NULL){ErrMsg(8.0);}
    
    if (No_reads == -1){
      if (readcov == -1.0){
        ErrMsg(2.0);
      }
    }

    int FixedSize = mypars->Length;
    const char* Sizefile = mypars->LengthFile;
    int meanlength;

    if (Sizefile == NULL){
      if (FixedSize == -1){ErrMsg(4.0);}
      else{
        meanlength = FixedSize;
      } 
    }

    if (Sizefile != NULL){
      double sum = 0; int n = 0;
      double* Length_tmp; double* Frequency_tmp;
      Length_tmp = new double[LENS];
      Frequency_tmp = new double[LENS];
      char buf[LENS];
      gzFile gz = Z_NULL;
      gz = gzopen(Sizefile,"r");
      assert(gz!=Z_NULL);
      while(gzgets(gz,buf,LENS)){
        Length_tmp[n] = atof(strtok(buf,"\n\t "));
        Frequency_tmp[n] = atof(strtok(NULL,"\n\t "));
        sum = sum + (Length_tmp[n]*(Frequency_tmp[n]-Frequency_tmp[n-1]));
        n++;
      }
      gzclose(gz);

      delete[] Frequency_tmp;delete[] Length_tmp;
      meanlength = (int) sum;

      if (FixedSize != -1){ErrMsg(5.0);}
    }

    faidx_t *seq_ref = NULL;
    seq_ref  = fai_load(fastafile);
    
    assert(seq_ref!=NULL);
    
    int chr_total = faidx_nseq(seq_ref);
    int Glob_seed = mypars->Glob_seed; 
    int threads1 = mypars->threads1;
    int threads2;
    if (mypars->threads2 == -1){threads2 = 1;}
    else{threads2 = mypars->threads2;}
    
    int Thread_specific_Read;
    if (No_reads != -1){
      if (No_reads == 1 && threads1 > 1){ErrMsg(2.6);}
     
      Thread_specific_Read = static_cast<int>(No_reads/threads1);
      if (readcov != -1){ErrMsg(3.0);}
    }
    else if (readcov != -1){
      size_t genome_size = 0;
      for (int i = 0; i < chr_total; i++){
        const char *chr_name = faidx_iseq(seq_ref,i);
        int chr_len = faidx_seq_len(seq_ref,chr_name);
        genome_size += chr_len;
      }
      Thread_specific_Read = static_cast<int>(((readcov*genome_size)/meanlength)/threads1);
    }

    fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,chr_total);
    fprintf(stderr,"\t-> Seed used: %d\n",Glob_seed);
    fprintf(stderr,"\t-> Number of threads used for sampling: %d and for writing down: %d\n",threads1,threads2);
    fprintf(stderr,"\t-> Number of simulated reads: %zd or coverage: %f\n",No_reads,readcov);

    const char* Adapt_flag;
    const char* Adapter_1;
    const char* Adapter_2;
    const char* Polynt;
    if (mypars->Adapter1 != NULL){
      Adapt_flag = "true";
      Adapter_1 = mypars->Adapter1;
      Adapter_2 = mypars->Adapter2;

      if (mypars->Poly != NULL){Polynt =mypars->Poly;}
      else{Polynt = "F";}
      
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

    const char* ErrorFlag;
    if (mypars->ErrorFlag != NULL){
      ErrorFlag = mypars->ErrorFlag;
      //fprintf(stderr,"ERROR FLAG IS %s\n",mypars->ErrorFlag);
    }
    else{ErrorFlag = "F";}

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
    }
    else{Briggs_Flag = "False";}
    
    const char* Specific_Chr[1024] = {};
    if (mypars->Chromosomes != NULL){
      int chr_idx_partial = 0;
      Specific_Chr[chr_idx_partial++] = strtok(strdup(mypars->Chromosomes),"\", \t");
      char *chrtok = NULL;
      while(((chrtok=strtok(NULL,"\", \t")))){
	      Specific_Chr[chr_idx_partial++] = strdup(chrtok);
	      assert(chr_idx_partial<MAXBINS);
      }
      Specific_Chr[chr_idx_partial++] = "\0";
    }
    
    //if(Specific_Chr[0]=='\0'){fprintf(stderr,"HURRA");}

    Create_se_threads(seq_ref,threads1,Glob_seed,Thread_specific_Read,filename,
                      Adapt_flag,Adapter_1,Adapter_2,OutputFormat,Seq_Type,
                      Param,Briggs_Flag,Sizefile,FixedSize,qualstringoffset,
                      QualProfile1,QualProfile2,threads2,QualStringFlag,Polynt,ErrorFlag,Specific_Chr);

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
  free((char *)mypars->Poly);
  //free((char *)mypars->ErrorFlag);
  delete mypars;

}
