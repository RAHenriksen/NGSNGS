#include <stdlib.h>
#include <random>
#include "mrand.h"

mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret = (mrand_t *) malloc(sizeof(mrand_t));
  ret->type = type_a;
#if defined(__linux__) || defined(__unix__)
  if(ret->type==0){
    fprintf(stderr,"In linux if -> drand48_data\n");
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);
    //i need to somehow print the value
  }
#elif defined(__APPLE__) || defined(__MACH__) //try to change this to linux and unix and then it works
  if(ret->type==1){
    fprintf(stderr,"In Apple loop if ->  APPLE LOOP\n");
    ret->eng = std::default_random_engine(seedval);
    ret->distr = std::uniform_real_distribution<double>(0, 1);
  }
#else
  else{
    fprintf(stderr,"type: %d is not defined, maybe problem with compiler macros\n",type_a);
    exit(0);
  }
#endif
  return ret;
}

double mrand_pop(mrand_t *mr){
  double res;
  #if defined(__linux__) || defined(__unix__)
    if(mr->type==0){
      drand48_r((struct drand48_data*)&mr->buf0,&res);
    }
  #elif defined(__APPLE__) || defined(__MACH__) 
    if(mr->type==1){
      res =  mr->distr(mr->eng);
    }
  #endif
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
