#ifndef NTSUBMODELS_H
#define NTSUBMODELS_H
#include "mrand.h"
#include <htslib/kstring.h>
#include <map>
#include "htslib/bgzf.h"
#include <zlib.h>

void ErrorSub(double randval,char seqchar[], int pos);

void MisMatchFileArray(double* freqval,const char* filename,int &mismatchcyclelength,int &elements);

int MisMatchFile(char seq[],mrand_t *mr,double* freqval,int LEN);

int MisMatchFile_kstring(kstring_t* seq,mrand_t *mr,double* freqval,int LEN);

typedef struct{
  int nreads;//this is nalignements
  double *fwD;
  double *bwD;
}mydataD;

std::map<int, mydataD> load_mismatch(const char* fname,int &printlength);

void parse_mismatch_Data(mydataD &md,double **dat,int howmany);

void MisMatchMetaFileArray(double* freqval,const char* filename,int &mismatchcyclelength,int &num_elem,const char* fileoutname);


#endif
