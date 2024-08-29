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
  int simmode;              //simulation mode
  int SamplThreads;           //sampling threads, used internally by pthread_create
  int CompressThreads;        //compression threads, used external by bgzf_mt and set_hts_options
  size_t nreads;              //Number of reads to simulate
  double coverage;            //Depth of coverage across the entire genome, estimated from rlen, flen, genomsize
  int Glob_seed;              //Local seeds are computed from the global. Only one seed needs to be supplied
  int Glob_seed_binary;       //Local seeds are computed from the global. Only one seed needs to be supplied
  outputformat_e OutFormat;  //fq, fq.gz, fa, fa.gz, sam, bam, cram
  char *OutName;              //prefix for output name
  char *VCFDumpFile;             //prefix for internal file recording potential variations to reference genome
  char *IndelDumpFile;        //prefix for internal file recording of sequencing errors specific for indels
  char *Reference;            //full filename for reference fasta
  seqtype_e seq_type;         //singleend or paired end.
  char *Adapter1;             //actual adapter sequence, R1, not flipped, reversed or completemented
  char *Adapter2;             //actual adapter sequence, R2, not flipped, reversed or completemented
  char *QualProfile1;         //filename for quality prof, for R1 or SE
  char *QualProfile2;         //filename for quality prof, for R2 used only in PE
  int FixedQual;
  char *SubProfile;           //filename for misincorperation, typespecific and position specific
  char *MisMatchMatrix_bdam;
  char *M3outname;
  int DoSeqErr;               //adding potential sequencing errors
  char *Briggs;               //the four briggs parameters for the none biotin model
  char *BriggsBiotin;         //the four briggs parameters in relation to Biotin
  int Duplicates;             //for the none biotin model, four potential fragments with independent deamination pattern are possible to generate, default PCR=1
  int CycleLength;            //cycle length with the maximum number of nucleotides to be generated, independent of the fragment length
  int LowerLimit;             //lower limit of fragment lengths, with some distributions perhaps having values in the pdf below 0, default of 30 due to the deamination
  int Length;                 //fragment length when fixed
  char *LengthFile;           //filename for distribution of frag lengths
  char *LengthDist;           //name of pdf used for simulating fraglengths (incl parameters)
  char *Poly;                 //poly-X tail, added after the adapter sequence if output reade are lower than inferred cycle length, e.g. -p G or -p A
  char *Chromosomes;          //subsetting chromosome of interested, -chr chrX
  char *BedFile;
  size_t flankingregion;
  int rng_type;               //pseudo-random number generator type, drand48 or rand or drand48_r etc
  char *vcffile;              //filename for bcf
  char *CommandRun;           //actual command run in same order
  int HeaderIndiv;            //samplename from VCF/BCF file
  char *NameIndiv;            //samplename from VCF/BCF file
  int Align;                  //Storing sequence reads with- or without alignment information in the sequence alignment map/format
  size_t KstrBuf;             //buffer size for kstring length
  char *Indel;                //adding stochastic indels
  double mutationrate;        //fixated mutation rate to alter reference genome
  int generations;           //mutation rate for a certain number of generations which alters the reference genome
  size_t referencevariations; //adding a fixed number of variations to the reference genome
  int MaskBed;
  int CaptureVCF;              //filename for bcf
  int linkage;
}argStruct;
argStruct *getpars(int argc,char ** argv);
void argStruct_destroy(argStruct *mypars);
#endif
