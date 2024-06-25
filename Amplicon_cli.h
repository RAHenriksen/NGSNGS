#ifndef AMPLICONCLI_H
#define AMPLICONCLI_H
#include <cstdlib>
#include <ctime>
#include <cstring>
#include "HelpPage.h"

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
}argStruct;
argStruct *amplicongetpars(int argc,char ** argv);

#endif
