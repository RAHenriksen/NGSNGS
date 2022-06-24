#include <stdlib.h>
#include <random>
#include <iostream>
#include "mrand.h"

//((double) rand_r(&seed)/ RAND_MAX);
mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret = (mrand_t *) malloc(sizeof(mrand_t));
  ret->type = type_a;

  if(ret->type==0){
    //fprintf(stderr,"In linux if -> drand48_data\n");
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);
    //i need to somehow print the value
  }
  if(ret->type==1){
    //fprintf(stderr,"In Apple loop if ->  APPLE LOOP\n");
    ret->eng = std::default_random_engine(seedval);
    ret->distr = std::uniform_real_distribution<float>(0, 1);
  }
  if(ret->type==2){
    //fprintf(stderr,"rand_r loop\n");
    ret->rand_r_seed = (unsigned int) seedval;
  }
  return ret;
}

double mrand_pop(mrand_t *mr){
  double res;
  if(mr->type==0){
    drand48_r((struct drand48_data*)&mr->buf0,&res);
  }
  else if(mr->type==1){
    res =  mr->distr(mr->eng);
  }
  else if(mr->type==2){
    res = (double) rand_r(&mr->rand_r_seed)/RAND_MAX;
  }
  else{
    fprintf(stderr,"Random parameter %d is not supported\n",mr->type);
    exit(0);
  }
  return res;
}


#ifdef __WITH_MAIN__

int main(){
  mrand_t *myrand;
  myrand = mrand_alloc(0,10);
  for(int i=0;i<10;i++)
    fprintf(stderr,"%d) : %f\n",i,mrand_pop(myrand));
  fprintf(stderr,"----------------\n");
  myrand = mrand_alloc(1,10);
  for(int i=0;i<10;i++)
    fprintf(stderr,"%d) : %f\n",i,mrand_pop(myrand));
  fprintf(stderr,"----------------\n");
  myrand = mrand_alloc(2,10);
  unsigned int seed = 10;
  for(int i=0;i<10;i++){
    //fprintf(stderr,"Rand %f \n",(double) rand_r(&seed)/RAND_MAX);
    fprintf(stderr,"%d) : %f\n",i,mrand_pop(myrand));
  }
  return 0;
}

#endif

//g++ mrand.cpp -std=c++11 -lm -lz -D__WITH_MAIN__ -o Rand
