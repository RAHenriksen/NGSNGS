#ifndef SAMPLE_QSCORES_H
#define SAMPLE_QSCORES_H
#include "mrand.h"
#include "RandSampling.h"
ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile);
void sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset);


#endif
