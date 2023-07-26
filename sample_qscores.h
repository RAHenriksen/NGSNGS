#ifndef SAMPLE_QSCORES_H
#define SAMPLE_QSCORES_H
#include "mrand.h"
#include "RandSampling.h"
ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile);
int sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset);
int sample_qscores_fix(char *bases,char *qscores, int qscorevalue,int len,mrand_t *mr,int simError, int ntcharoffset);

#endif
