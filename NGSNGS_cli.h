#ifndef NGSNGSCLI_H
#define NGSNGSCLI_H
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "HelpPage.h"
enum outputformat_e {unknownT, faT, fagzT, fqT, fqgzT, samT, bamT,cramT };
enum seqtype_e {unknownTT, SE,PE};
enum dist_e {unknownDist,Uni,Norm,LogNorm,Pois,Exp,Gam};
enum poly_nt {unknownNt,A,G,C,T,N};

typedef struct{
  int SamplThreads; //sampling threads, used internally by pthread_create
  int CompressThreads; //compression threads, used external by bgzf_mt and set_hts_options
  size_t nreads; // Thorfinn made mistake
  double coverage; // Coverage of breadth, estimated from rlen, flen, genomsize, superfancy
  int Glob_seed; // local seeds are computed determistly from global. Only one seed needs to be supplied
  outputformat_e OutFormat ; //fq, fq.gz, fa, fa.gz, sam, bam, cram, ultrafancy
  char *OutName; //prefix for output, final will be determined by [OutName.OutFormat]
  char *Reference; //full filename for reference fasta
  seqtype_e seq_type; // singleend or paired end.
  char *Adapter1; //actual adapterseq, R1, not flipped, reversed or completemented
  char *Adapter2; //actual adapterseq, R2, not flipped, reversed or completemented
  char *QualProfile1;// filename for quality prof, for R1 or SE
  char *QualProfile2;// filename for quality prof, for R2 used only in PE
  char *SubProfile;// filename for misincorperation, typespecific and position specific
  int DoSeqErr;// should we not add errors?, {T,F}
  char *Briggs; // the four briggs parameters
  int Length; // fragment length when fixed
  char *LengthFile; //filename for distribution of frag lengths
  char *LengthDist;//name of pdf used for simulating fraglengths (incl parameters)
  char *Poly; // should poly be added possible values -p G or -p A
  char *Chromosomes; //for subsetting chromosome of interested, -chr chrX
  int rng_type; //RNG type, drand48 or rand or drand48_r etc
  char *vcffile; // filename for bcf
  char *CommandRun; // actual command run in same order
  char* HeaderIndiv; //samplename from VCF/BCF file
  int NoAlign;// This option is cool, but explaining it takes up to much space in comment
  size_t KstrBuf; // Buffer size for kstring length
}argStruct;
argStruct *getpars(int argc,char ** argv);
void argStruct_destroy(argStruct *mypars);
#endif
