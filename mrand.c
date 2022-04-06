#include <stdlib.h>
#include "mrand.h"

mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret = (mrand_t *) malloc(sizeof(mrand_t));
  ret->type = type_a;
  if(ret->type==0)
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);

  return ret;
}

double mrand_pop(mrand_t *mr){
  double res;
  if(mr->type==0)
    drand48_r((struct drand48_data*)&mr->buf0,&res);

  return res;
}
