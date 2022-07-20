#include <stdlib.h>
#include <random>
#include <iostream>
#include "mrand.h"

mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret = (mrand_t *) malloc(sizeof(mrand_t));
  ret->type = type_a;

#ifdef __APPLE__
  if(ret->type==0){
    fprintf(stderr,"\t-> Problem with drand48 reentrant, will default to erand48\n");
    ret->type = 3;
  }
#else
  if(ret->type==0)
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);
#endif
  if(ret->type==1){
    ret->eng = std::default_random_engine(seedval);
    ret->distr = std::uniform_real_distribution<float>(0, 1);
  }
  if(ret->type==2)
    ret->rand_r_seed = (unsigned int) seedval;
  if(ret->type==3){
    ret->rand_r_seed = (unsigned int) seedval;
    ret->xsubi[2] = ret->rand_r_seed >> 16;
    ret->xsubi[1] = ret->rand_r_seed & 0xffffl;
    ret->xsubi[0] = 0x330e;
  }
  return ret;
}

double mrand_pop(mrand_t *mr){
  double res;
  if(mr->type==0){
    #if defined(__linux__) || defined(__unix__)
    drand48_r((struct drand48_data*)&mr->buf0,&res);
    #endif
  }
  else if(mr->type==1){
    res =  mr->distr(mr->eng);
  }
  else if(mr->type==2){
    res = (double) rand_r(&mr->rand_r_seed)/RAND_MAX;
  }
  else if(mr->type==3){
    res = erand48(mr->xsubi);
  }
  else{
    fprintf(stderr,"Random parameter %d is not supported\n",mr->type);
    exit(0);
  }
  return res;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  mrand_t *myrand;
  int type =0;
  int seed =1;
  int nitems =10;
  if(argc==4){
    type = atoi(argv[1]);
    seed = atoi(argv[2]);
    nitems = atoi(argv[3]);
  }
  fprintf(stderr,"type: %d seed: %d nitems: %d\n",type,seed,nitems);
  myrand = mrand_alloc(type,seed);
  for(int i=0;i<nitems;i++)
    fprintf(stdout,"%d) : %f\n",i,mrand_pop(myrand));
  return 0;
}

#endif

//g++ mrand.cpp -std=c++11 -lm -lz -D__WITH_MAIN__ -o Rand
