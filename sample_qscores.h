#ifndef SAMPLE_QSCORES_H
#define SAMPLE_QSCORES_H
#include "mrand.h"
#include "RandSampling.h"
#include <htslib/kstring.h>

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,int &inferredreadcycle);
int sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset);
int sample_qscores_amplicon(kstring_t* seq,kstring_t* qual,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset);

int sample_qscores_fix(char *bases,char *qscores, int qscorevalue,int len,mrand_t *mr,int simError, int ntcharoffset);
int sample_qscores_fix_amplicon(mrand_t *mr,kstring_t* seq,int qscoresval, int ntcharoffset);

#endif
