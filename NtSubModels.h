#ifndef NTSUBMODELS_H
#define NTSUBMODELS_H
#include "mrand.h"

void ErrorSub(double randval,char seqchar[], int pos);
double* MisMatchFileArray(double* freqval,const char* filename,int &mismatchcyclelength);
int MisMatchFile(char seq[],mrand_t *mr,double* freqval,int LEN);

#endif
