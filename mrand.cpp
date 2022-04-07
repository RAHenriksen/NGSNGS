#include <stdlib.h>
#include <random>
#include "mrand.h"

mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret = (mrand_t *) malloc(sizeof(mrand_t));
  ret->type = type_a;
#ifndef __APPLE__
  if(ret->type==0)
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);
  else
#endif
    if(ret->type==1){
    ret->eng = std::default_random_engine(seedval);
    ret->distr = std::uniform_real_distribution<float>(0, 1);

  } else{
    fprintf(stderr,"type: %d is not defined, maybe problem with compiler macros\n",type_a);
    exit(0);
  }
  return ret;
}

double mrand_pop(mrand_t *mr){
  double res;
  if(mr->type==0)
    drand48_r((struct drand48_data*)&mr->buf0,&res);
  if(mr->type==1){
    res =  mr->distr(mr->eng);

  }
  return res;
}


#ifdef __WITH_MAIN__

int main(){
  mrand_t *myrand = mrand_alloc(1,10);
  for(int i=0;i<10;i++)
    fprintf(stderr,"%d) : %f\n",i,mrand_pop(myrand));
  return 0;
}

#endif
