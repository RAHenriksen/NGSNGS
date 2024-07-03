#ifndef AMPLICONCLI_H
#define AMPLICONCLI_H
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "HelpPage.h"

enum ampliconformat_e {unknownT, faT, fagzT, fqT, fqgzT, samT, bamT};

// Define parameters
typedef struct{
  const char* Amplicon_in_pars;
  const char* Amplicon_out_pars;
  char *BriggsBiotin;  //the four briggs parameters in relation to Biotin
  char *Indel;         //the four indel parameters in relation to random indels
  int Threads;
  long int Seed;
  char *SubProfile;           //filename for misincorperation, typespecific and position specific
  int rng_type;
  ampliconformat_e OutFormat ;  //fq, fq.gz, fa, fa.gz, sam, bam, cram
  char *CommandRun;           //actual command run in same order

}argStruct;

argStruct *amplicongetpars(int argc,char ** argv);

void amplicongetpars_destroy(argStruct *mypars);

#endif
